
__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"

import statistics

import pysam


def read_cnvkit_cns(input_cns):
    """Read a CNVkit .cns file into {chrom: [[start, end, log2], ...]}, sorted by start position."""
    cns_dict = {}
    with open(input_cns) as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_i = header.index("chromosome")
        start_i = header.index("start")
        end_i = header.index("end")
        log2_i = header.index("log2")
        for line in f:
            if not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            chrom = columns[chrom_i]
            start = int(columns[start_i])
            end = int(columns[end_i])
            log2 = float(columns[log2_i])
            cns_dict.setdefault(chrom, []).append([start, end, log2])
    for chrom in cns_dict:
        cns_dict[chrom].sort(key=lambda seg: seg[0])
    return cns_dict


def lookup_local_log2(cns_dict, chrom, pos):
    """Return the log2 ratio of the segment covering pos, or None if uncovered."""
    for start, end, log2 in cns_dict.get(chrom, []):
        if start <= pos <= end:
            return log2
    return None


def infer_sex_from_cnr(input_cnr):
    """
    Infer sample sex from a CNVkit .cnr file by comparing this sample's own
    median chrX log2 to its own median autosomal log2. Defaults to "female"
    (chrX treated as diploid) when there's too little chrX data to tell.
    """
    autosome_log2 = []
    chrx_log2 = []
    with open(input_cnr) as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_i = header.index("chromosome")
        log2_i = header.index("log2")
        for line in f:
            if not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            chrom = columns[chrom_i].lstrip("chr")
            log2 = float(columns[log2_i])
            if chrom == "X":
                chrx_log2.append(log2)
            elif chrom not in ("Y", "MT", "M"):
                autosome_log2.append(log2)

    if len(chrx_log2) < 10 or not autosome_log2:
        return "female"

    diff = statistics.median(chrx_log2) - statistics.median(autosome_log2)
    return "male" if diff < -0.5 else "female"


def read_cns_segments(cns_path):
    """
    Read a CNVkit purity/BAF-aware .loh.cns file into {chrom: [[start, end, baf_or_None], ...]}.
    baf is None where CNVkit had too few germline SNPs to fit one for that segment
    (routinely empty on a targeted panel, not a rare edge case).
    """
    cns_dict = {}
    with open(cns_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_i = header.index("chromosome")
        start_i = header.index("start")
        end_i = header.index("end")
        baf_i = header.index("baf")
        for line in f:
            if not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            chrom = columns[chrom_i]
            start = int(columns[start_i])
            end = int(columns[end_i])
            baf = float(columns[baf_i]) if columns[baf_i] != "" else None
            cns_dict.setdefault(chrom, []).append([start, end, baf])
    for chrom in cns_dict:
        cns_dict[chrom].sort(key=lambda seg: seg[0])
    return cns_dict


def lookup_local_baf(cns_dict, chrom, pos):
    """Return the fitted baf of the segment covering pos, or None if uncovered/unfit."""
    for start, end, baf in cns_dict.get(chrom, []):
        if start <= pos <= end:
            return baf
    return None


def is_likely_germline(af, baf, af_germline_lower_limit, af_germline_upper_limit, baf_tolerance):
    """
    Without local BAF, use the fixed diploid-heterozygous window (~0.5). When
    CNVkit fitted a reliable baf for the covering segment, a germline SNP
    there clusters near baf or its mirror (1-baf) instead of near 0.5 - e.g.
    a segment under LOH/allelic imbalance shifts a germline het's AF well
    away from 50% while still being germline, which the fixed window alone
    would miss.
    """
    if baf is None:
        return af_germline_lower_limit <= af <= af_germline_upper_limit
    return abs(af - baf) <= baf_tolerance or abs(af - (1 - baf)) <= baf_tolerance


# Below this purity, local CN reads aren't trustworthy (lab convention:
# ~10-20% for FFPE) - also used as the bisection search floor below.
PURITY_FLOOR = 0.10
# Below this |log2ratio|, the signal is indistinguishable from ordinary
# CNVkit segment-level noise and is treated as neutral regardless of purity.
LOG2_NOISE_FLOOR = 0.1
BISECTION_TOLERANCE = 1e-6
BISECTION_MAX_ITERATIONS = 60


def local_cn_at_purity(log2ratio, normal_cn, purity, purity_floor=PURITY_FLOOR):
    """
    Invert the log2 dilution equation for tumor copy number at a given
    purity: 2^log2ratio * normal_cn = purity*CN_t + (1-purity)*normal_cn.

    Floored at 1 (not 0): a variant was actually called here, so at least
    one tumor-derived copy must exist - without this, cn_t->0 flips the
    correction upward as a model artifact rather than a real signal.

    Kept continuous, not rounded - see correct_vaf_for_copy_number.
    """
    p = max(purity, purity_floor)
    cn_t = normal_cn * (2 ** log2ratio - (1 - p)) / p
    return max(cn_t, 1.0)


def _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, purity, purity_floor=PURITY_FLOOR):
    """
    Given a candidate purity, de-dilute log2ratio into local copy number,
    derive the assumed mutant-copy count, and solve the VAF equation for
    the tumor fraction that implies.

    Mutant allele is assumed amplified on a gain (m = cn_t - (normal_cn-1))
    and to be the retained copy on a loss/LOH (m = cn_t); reduces to the
    plain heterozygous model (m=1) at neutral CN.

    return: (cn_t, m, adjusted_tc)
    """
    if log2ratio is None:
        cn_t = float(normal_cn)
    else:
        cn_t = local_cn_at_purity(log2ratio, normal_cn, purity, purity_floor)

    # Tolerance avoids floating-point noise (e.g. 1.9999999999996 instead of
    # 2.0) flipping this into the loss branch right at the neutral boundary.
    if cn_t >= normal_cn - 1e-9:
        m = cn_t - (normal_cn - 1)
    else:
        m = cn_t
    m = max(m, 1e-9)  # avoid division by zero on a homozygous deletion (cn_t=0)

    denominator = m - vaf * (cn_t - normal_cn)
    if denominator <= 0:
        adjusted_tc = 1.0  # observed VAF implies >=100% purity under this model
    else:
        adjusted_tc = (vaf * normal_cn) / denominator
        adjusted_tc = min(max(adjusted_tc, 0.0), 1.0)

    return cn_t, m, adjusted_tc


def correct_vaf_for_copy_number(
    vaf, log2ratio, normal_cn=2, log2_noise_floor=LOG2_NOISE_FLOOR, purity_floor=PURITY_FLOOR
):
    """
    Convert a somatic VAF into a tumor-fraction estimate corrected for local
    copy number (vs. today's plain VAF*2). See _cn_and_adjusted_tc for the
    allele-count model.

    cn_t depends on purity, which is what we're estimating - two equations
    (log2ratio, vaf), one unknown (purity), solved by bisection on
    h(p) = adjusted_tc(cn_t(p)) - p. cn_t is kept continuous throughout the
    solve and rounded only at the end, for reporting: rounding it mid-solve
    turns h into a step function with spurious extra roots near rounding
    boundaries (confirmed empirically).

    Falls back to the uncorrected estimate when no self-consistent purity
    exists (h single-signed across [purity_floor, 1] - typically an
    aggressive gain incompatible with the VAF at any purity). Does not
    correct for subclonality (VAF reflects purity*CCF, not separable from
    purity alone given a single variant).

    Below purity_floor, the correction is evaluated assuming exactly
    purity_floor rather than skipped or searched for, and reported
    unfloored (neutral CN reduces to raw_tc unchanged; a real gain/loss
    still gets a genuine correction, which can land either side of the
    floor).

    return: (raw_tc, adjusted_tc, cn_t, mutant_copies)
            cn_t is a rounded label only - adjusted_tc/mutant_copies come
            from the continuous solve, not from cn_t
    """
    raw_tc = vaf * 2

    if log2ratio is None or abs(log2ratio) < log2_noise_floor:
        cn_t, m, adjusted_tc = _cn_and_adjusted_tc(vaf, None, normal_cn, 1.0)
        return raw_tc, adjusted_tc, round(cn_t), m

    if raw_tc < purity_floor:
        # Too low to search for a self-consistent purity - evaluate the
        # correction assuming exactly purity_floor instead, unfloored.
        cn_t, m, adjusted_tc = _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, purity_floor, purity_floor)
        return raw_tc, adjusted_tc, round(cn_t), m

    def h(p):
        return _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, p, purity_floor)[2] - p

    lo, hi = purity_floor, 1.0
    h_lo, h_hi = h(lo), h(hi)

    if abs(h_lo) < BISECTION_TOLERANCE:
        p_star = lo
    elif abs(h_hi) < BISECTION_TOLERANCE:
        p_star = hi
    elif (h_lo > 0) == (h_hi > 0):
        # No sign change - no self-consistent purity, fall back to raw_tc
        return raw_tc, raw_tc, normal_cn, 1.0
    else:
        p_star = (lo + hi) / 2
        for _ in range(BISECTION_MAX_ITERATIONS):
            p_star = (lo + hi) / 2
            h_mid = h(p_star)
            if abs(h_mid) < BISECTION_TOLERANCE or (hi - lo) < BISECTION_TOLERANCE:
                break
            if (h_mid > 0) == (h_lo > 0):
                lo, h_lo = p_star, h_mid
            else:
                hi, h_hi = p_star, h_mid

    # adjusted_tc/m come from the continuous solve, never re-derived from
    # the rounded cn_t below (that would reintroduce the same instability).
    cn_t_continuous, m, adjusted_tc = _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, p_star, purity_floor)
    cn_t = round(cn_t_continuous)

    return raw_tc, adjusted_tc, cn_t, m


def read_snv_vcf_and_find_max_af(input_snv_vcf, filter_dict):
    snv_vcf = pysam.VariantFile(input_snv_vcf)

    best_variant = []

    vep_fields = {}
    for record in snv_vcf.header.records:
        if record.type == "INFO":
            if record['ID'] == "CSQ":
                vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split('">')[0].split("|"))}

    for record in snv_vcf.fetch():
        if "COMPLEXAF" in record.info:
            # VarDict's synthetic pseudo-record for one component of a decomposed
            # complex variant, not an independently-supported call - it carries no
            # QUAL and none of the usual per-record QC annotation (NM/PMEAN/SN/...),
            # so it can slip through filter_dict undetected rather than genuinely
            # passing QC. The real, fully-annotated call for the same event is a
            # separate record at the same position.
            continue

        if "AF" not in record.info:
            # A synthetic codon-level substitution record (INFO only ever
            # carries AA/Artifact/CSQ) from this VCF's codon_snvs merge step -
            # not an independently-called variant, so it has no AF/DP/QC
            # annotation of its own to filter or estimate purity from. Without
            # this guard, only ever surviving to the AF lookup below by luck
            # of which filters happen to be configured (e.g. currently only
            # the Artifact=-1 sentinel catches it) crashes with a KeyError.
            continue

        vep = record.info["CSQ"][0]
        vep_dict = dict(zip(vep_fields.keys(), vep.split("|")))

        filtered = False
        for filter in filter_dict:
            if filter in record.info:
                if filter == "Artifact":
                    a1 = int(record.info[filter][0])
                    a2 = int(record.info[filter][1])
                    if a1 > filter_dict[filter][1] or a2 > filter_dict[filter][1] or a1 == -1 or a2 == -1:
                        filtered = True
                elif filter == "AF":
                    if record.info[filter][0] > filter_dict[filter][1]:
                        filtered = True
                elif filter_dict[filter][0] == "min":
                    if record.info[filter] < filter_dict[filter][1]:
                        filtered = True
                elif filter_dict[filter][0] == "max":
                    if record.info[filter] > filter_dict[filter][1]:
                        filtered = True
                elif filter_dict[filter][0] == "present":
                    if filter_dict[filter][1] not in record.info[filter]:
                        filtered = True
            elif filter in vep_dict:
                if vep_dict[filter] == "":
                    continue
                if filter_dict[filter][0] == "min":
                    if float(vep_dict[filter]) < filter_dict[filter][1]:
                        filtered = True
                elif filter_dict[filter][0] == "max":
                    if float(vep_dict[filter]) > filter_dict[filter][1]:
                        filtered = True
                elif filter_dict[filter][0] == "exact":
                    if vep_dict[filter] in filter_dict[filter][1]:
                        filtered = True
            elif filter == "Other":
                if not (vep_dict["IMPACT"] == "HIGH" or
                        vep_dict["Existing_variation"].count("COSV") > 1 or
                        vep_dict["CLIN_SIG"].find("drug_response") != -1 or
                        vep_dict["CLIN_SIG"].find("pathogenic") != -1 or
                        ("Hotspot" in record.info and record.info["Hotspot"] == "1-hotspot")
                        ):
                    filtered = True
            elif filter == "CHIP_genes":
                if vep_dict["SYMBOL"] in filter_dict[filter][1]:
                    filtered = True

        if not filtered:
            best_variant.append([record.info["AF"][0], record.chrom, record.pos, str(record)])

    best_variant.sort(key=lambda x: x[0], reverse=True)
    return best_variant


def drop_likely_germline(candidates, loh_cns_dict, af_germline_lower_limit, af_germline_upper_limit, baf_tolerance):
    """Drop candidates whose AF matches the local CN-aware germline pattern - see is_likely_germline."""
    kept = []
    for af, chrom, pos, record_str in candidates:
        baf = lookup_local_baf(loh_cns_dict, chrom, pos)
        if not is_likely_germline(af, baf, af_germline_lower_limit, af_germline_upper_limit, baf_tolerance):
            kept.append([af, chrom, pos, record_str])
    return kept


def write_tc(output_tc, raw_tc, adjusted_tc, raw_tc_all=None, adjusted_tc_all=None):
    '''
    Write the raw/adjusted driver-based TC, plus the passenger-based
    raw_tc_all/adjusted_tc_all ("NA" if not promoted - see __main__).
    Returns the output string to simplify unit testing.
    '''
    output = open(output_tc, "w")
    output.write(
        "Percentage ctDNA based on driver SNVs (raw)\tPercentage ctDNA based on driver SNVs (adjusted)\t"
        "Percentage ctDNA based on all SNVs (raw)\tPercentage ctDNA based on all SNVs (adjusted)\n"
    )
    raw_tc_all_str = f"{raw_tc_all*100:.1f}%" if raw_tc_all is not None else "NA"
    adjusted_tc_all_str = f"{adjusted_tc_all*100:.1f}%" if adjusted_tc_all is not None else "NA"
    line = f"{raw_tc*100:.1f}%\t{adjusted_tc*100:.1f}%\t{raw_tc_all_str}\t{adjusted_tc_all_str}\n"
    output.write(line)
    output.close()
    return line


def write_ctDNA_fraction_info(output_file_name, snv_info_list):
    '''
    Write the SNV candidates used to estimate TC, with their copy-number-
    correction details.

    param snv_info_list: [[raw_tc, adjusted_tc, cn_t, mutant_copies, normal_cn, source, VCF_record], ...]
                          source is "driver" or "passenger"
    '''
    output = open(output_file_name, "w")
    output.write("SNVs passing all filtering\n")
    output.write("raw_%\tadjusted_%\tlocal_CN_t\tassumed_mutant_copies\tnormal_CN_used\tsource\tVCF_record\n")
    for raw_tc, adjusted_tc, cn_t, m, normal_cn, source, record_str in snv_info_list:
        output.write(f"{raw_tc*100:.1f}%\t{adjusted_tc*100:.1f}%\t{cn_t}\t{m:.2f}\t{normal_cn}\t{source}\t{record_str}")
    output.close()


if __name__ == "__main__":
    input_vcf = snakemake.input.vcf
    input_cnvkit_cns = snakemake.input.cnvkit_cns
    input_cnvkit_cnr = snakemake.input.cnvkit_cnr
    input_loh_cns = snakemake.input.loh_cns
    output_ctDNA_fraction = snakemake.output.ctDNA_fraction
    output_ctDNA_fraction_info = snakemake.output.ctDNA_fraction_info

    callers = snakemake.params.callers
    if isinstance(callers, list):
        callers = callers[0]

    filter_dict = {
        "PositionNrSD": ["min", snakemake.params.min_position_nr_sd],
        "PanelMedian": ["max", snakemake.params.max_panel_median],
        "Artifact": ["max", snakemake.params.artifact_limit],
        "CALLERS": ["present", callers],
        "MQ": ["min", snakemake.params.min_mq],
        "MSI": ["max", snakemake.params.max_msi],
        "NM": ["max", snakemake.params.max_nm],
        "ODDRATIO": ["max", snakemake.params.max_odd_ratio],
        "PMEAN": ["min", snakemake.params.min_pmean],
        "QUAL": ["min", snakemake.params.min_qual],
        "SBF": ["min", snakemake.params.min_sbf],
        "AF": ["max", snakemake.params.max_af],
        "MAX_AF": ["max", snakemake.params.max_gnomad_af],
        "Consequence": ["exact", snakemake.params.excluded_consequences],
        "CHIP_genes": ["", snakemake.params.chip_genes],
        "Other": ["", []]
    }

    driver_candidates = read_snv_vcf_and_find_max_af(input_vcf, filter_dict)

    # Same filters minus the "Other" driver gate and the Consequence gate, so
    # passenger (incl. synonymous/UTR) mutations are eligible too - a less
    # trustworthy but still useful fallback/cross-check; QC filters still apply.
    filter_dict_all = {name: spec for name, spec in filter_dict.items() if name not in ("Other", "Consequence")}
    all_candidates = read_snv_vcf_and_find_max_af(input_vcf, filter_dict_all)

    # Drop germline calls that only the fixed AF<max_af window would miss -
    # e.g. a known benign SNP whose AF is shifted well away from 50% by local
    # LOH/allelic imbalance, which dropping the Other/Consequence gates above
    # would otherwise let through as a "passenger".
    loh_cns_dict = read_cns_segments(input_loh_cns)
    driver_candidates = drop_likely_germline(
        driver_candidates, loh_cns_dict, snakemake.params.af_germline_lower_limit,
        snakemake.params.af_germline_upper_limit, snakemake.params.cn_baf_tolerance
    )
    all_candidates = drop_likely_germline(
        all_candidates, loh_cns_dict, snakemake.params.af_germline_lower_limit,
        snakemake.params.af_germline_upper_limit, snakemake.params.cn_baf_tolerance
    )

    cns_dict = read_cnvkit_cns(input_cnvkit_cns)
    inferred_sex = infer_sex_from_cnr(input_cnvkit_cnr)

    def correct_candidate(af, chrom, pos):
        normal_cn = 1 if (chrom.lstrip("chr") == "X" and inferred_sex == "male") else 2
        log2ratio = lookup_local_log2(cns_dict, chrom, pos)
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(
            af, log2ratio, normal_cn, snakemake.params.log2_noise_floor, snakemake.params.purity_floor
        )
        return raw_tc, adjusted_tc, cn_t, m, normal_cn

    snv_info_list = []
    for af, chrom, pos, record_str in driver_candidates:
        raw_tc, adjusted_tc, cn_t, m, normal_cn = correct_candidate(af, chrom, pos)
        snv_info_list.append([raw_tc, adjusted_tc, cn_t, m, normal_cn, "driver", record_str])

    if snv_info_list:
        raw_tc, adjusted_tc = snv_info_list[0][0], snv_info_list[0][1]
    else:
        raw_tc, adjusted_tc = 0, 0

    # Only the best passenger candidate is evaluated; surfaced as
    # raw_tc_all/adjusted_tc_all only when no driver was found or it exceeds
    # the driver's raw_tc.
    raw_tc_all, adjusted_tc_all = None, None
    if all_candidates:
        af, chrom, pos, record_str = all_candidates[0]
        candidate_raw_tc, candidate_adjusted_tc, cn_t, m, normal_cn = correct_candidate(af, chrom, pos)
        snv_info_list.append([candidate_raw_tc, candidate_adjusted_tc, cn_t, m, normal_cn, "passenger", record_str])
        if not driver_candidates or candidate_raw_tc > raw_tc:
            raw_tc_all, adjusted_tc_all = candidate_raw_tc, candidate_adjusted_tc

    write_tc(output_ctDNA_fraction, raw_tc, adjusted_tc, raw_tc_all, adjusted_tc_all)
    write_ctDNA_fraction_info(output_ctDNA_fraction_info, snv_info_list)
