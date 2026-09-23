
__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"

import statistics

import pysam


def read_cnvkit_cns(input_cns):
    """
    Read a CNVkit .cns segment file into a per-chromosome sorted list of
    [start, end, log2] for coordinate-based lookup.

    param input_cns: path to a CNVkit .cns file
    return: dict {chrom: [[start, end, log2], ...]}, sorted by start position
    """
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
    """
    Find the log2 ratio of the CNVkit segment covering a genomic position.

    param cns_dict: dict created by read_cnvkit_cns
    param chrom: chromosome name (must match the .cns file's naming)
    param pos: 1-based genomic position
    return: log2 ratio (float), or None if no segment covers this position
    """
    for start, end, log2 in cns_dict.get(chrom, []):
        if start <= pos <= end:
            return log2
    return None


def infer_sex_from_cnr(input_cnr):
    """
    Infer sample sex from a CNVkit .cnr bin-level file, by comparing this
    sample's own median chrX log2 to its own median autosomal log2
    (self-normalized - robust regardless of whether CNVkit's own reference
    already applies its own sex-aware chrX normalization).

    param input_cnr: path to a CNVkit .cnr file
    return: "male" if chrX reads ~1 copy relative to autosomes, else "female".
            "female" is also the fallback when there isn't enough chrX data to
            tell - the safe, non-destructive default, since it just means
            chrX gets treated like an autosome (normal_CN=2).
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


# Below this purity, CNVs aren't reliably callable in the first place (lab
# convention: ~20% for FFPE, though down to ~10% is considered usable) - so
# copy-number correction is skipped entirely below it, rather than trusting
# a local CN read at a purity where that read itself can't be trusted. Also
# used as the bisection search floor whenever correction does run, so the
# solver can never wander back down into that untrusted range either.
PURITY_FLOOR = 0.10
# Below this |log2ratio|, the reading is indistinguishable from ordinary
# CNVkit segment-level measurement noise (typically ~0.1-0.2) and is treated
# as neutral regardless of purity - a real signal is required to correct on,
# not just a low-enough purity to search at. Complementary to PURITY_FLOOR,
# not a replacement for it: PURITY_FLOOR gates on whether the *sample* is
# pure enough to trust any CN read; this gates on whether *this particular*
# log2ratio reading is itself distinguishable from noise. At low purity the
# CN-rounding boundary (see local_cn_at_purity) sits far inside this noise
# band, so without this floor, ordinary noise can look like a multi-copy
# gain and pull the corrected estimate down by 50%+ from noise alone.
LOG2_NOISE_FLOOR = 0.1
BISECTION_TOLERANCE = 1e-6
BISECTION_MAX_ITERATIONS = 60


def local_cn_at_purity(log2ratio, normal_cn, purity, purity_floor=PURITY_FLOOR):
    """
    Invert the log2 dilution equation for the tumor's own local copy number
    at a given purity:

        2^log2ratio * normal_cn = purity * CN_tumor + (1 - purity) * normal_cn

    At purity=1 this reduces to the purity-naive normal_cn * 2^log2ratio
    conversion. At lower purity, the same observed log2ratio implies a more
    extreme CN_tumor than the naive conversion reports, since the signal is
    diluted by normal-cell admixture. purity is floored at purity_floor
    before being used as a divisor, since it approaches 0 as purity -> 0.

    Floored at 1, not 0: this is only ever evaluated at a locus where a
    variant was actually called, so at least one tumor-derived copy must be
    present there - a true homozygous deletion (cn_t=0) couldn't produce any
    tumor-derived reads to call a variant from in the first place. Without
    this floor, a candidate cn_t approaching 0 drives the mutant-copy count m
    toward 0 too, which can flip the corrected purity estimate upward
    (sometimes sharply) instead of the usual downward correction - a model
    artifact, not a real biological signal, for a variant we know exists.

    Deliberately NOT rounded to an integer - see correct_vaf_for_copy_number
    for why the self-consistency solve needs this to stay continuous.

    param log2ratio: local CNVkit log2 copy ratio
    param normal_cn: copy number of this locus in a normal cell
    param purity: current purity estimate used to de-dilute log2ratio
    param purity_floor: see PURITY_FLOOR
    return: local tumor copy number (continuous, not rounded), >= 1
    """
    p = max(purity, purity_floor)
    cn_t = normal_cn * (2 ** log2ratio - (1 - p)) / p
    return max(cn_t, 1.0)


def _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, purity, purity_floor=PURITY_FLOOR):
    """
    Given a candidate purity, de-dilute log2ratio into a (continuous) local
    copy number (or fall back to the neutral model if log2ratio is None),
    derive the assumed mutant-copy count for that copy-number state, and
    solve the VAF equation for the tumor fraction that state implies.

    Model: the mutant allele is assumed to be the one preferentially
    amplified on a gain (m = CN_t - (normal_cn - 1), i.e. exactly
    normal_cn - 1 copies stay wild-type and everything else is mutant), and
    the retained copy/copies are assumed mutant on a loss/LOH (m = CN_t).
    At CN_t == normal_cn this reduces to the plain heterozygous model (m=1).

    return: (cn_t, m, adjusted_tc)
    """
    if log2ratio is None:
        cn_t = float(normal_cn)
    else:
        cn_t = local_cn_at_purity(log2ratio, normal_cn, purity, purity_floor)

    # Tolerance, not a plain >=: local_cn_at_purity's floating-point chain can
    # land a fraction below normal_cn at exact neutrality (e.g. 1.9999999999996
    # instead of 2.0), which would otherwise flip this into the loss branch
    # and silently produce the wrong m right at the neutral boundary.
    if cn_t >= normal_cn - 1e-9:
        m = cn_t - (normal_cn - 1)
    else:
        m = cn_t
    m = max(m, 1e-9)  # guard against division by zero on a homozygous deletion (cn_t=0)

    denominator = m - vaf * (cn_t - normal_cn)
    if denominator <= 0:
        # Degenerate case - the observed VAF implies >=100% purity under this model
        adjusted_tc = 1.0
    else:
        adjusted_tc = (vaf * normal_cn) / denominator
        adjusted_tc = min(max(adjusted_tc, 0.0), 1.0)

    return cn_t, m, adjusted_tc


def correct_vaf_for_copy_number(
    vaf, log2ratio, normal_cn=2, log2_noise_floor=LOG2_NOISE_FLOOR, purity_floor=PURITY_FLOOR
):
    """
    Convert a somatic VAF into a tumor-fraction estimate, correcting for
    local copy number instead of always assuming a plain diploid
    heterozygous locus (today's VAF*2 behaviour). See _cn_and_adjusted_tc
    for the allele-count model.

    CN_t itself depends on purity (see local_cn_at_purity), which is exactly
    what this function is trying to estimate. log2ratio and vaf are both
    fixed observations with a single unknown, the true purity p - two
    equations, one unknown - so this solves for the self-consistent p
    directly via bisection on h(p) = adjusted_tc(cn_t(p)) - p over
    p in [PURITY_FLOOR, 1], with cn_t kept CONTINUOUS (not rounded to an
    integer) throughout the solve. cn_t is rounded to the nearest integer
    only afterwards, purely for the returned/reported value.

    Rounding cn_t to an integer *during* the solve was tried and rejected:
    it turns h into a step function, and multiple "self-consistent" points
    can appear near rounding boundaries that have nothing to do with a real
    solution (confirmed empirically - at low purity, several adjacent
    integer CN states can each look self-consistent for the same observed
    data, with no principled way to pick the "right" one from a single
    variant's numbers alone). The continuous solve avoids this: the
    underlying problem is well-posed (2 equations, 1 unknown) and h is
    smooth, so plain bisection reliably finds the single correct root.

    Not every (vaf, log2ratio) pair has a self-consistent solution under this
    model - h(p) can be single-signed across the whole range (most often for
    an aggressive gain, where the "mutant preferentially amplified" copy
    assumption is simply incompatible with the observed VAF at any purity).
    When that happens there's no principled p to report, so this falls back
    to the uncorrected neutral model (adjusted_tc = raw_tc) rather than
    extrapolating to a boundary that could be wildly wrong either way.

    This does not correct for subclonality (a subclonal variant's VAF
    reflects purity * CCF, which this model can't separate from purity alone
    given only a single variant's VAF) - only for the copy-number dilution
    bias.

    Below PURITY_FLOOR, raw_tc is too low to search for a self-consistent
    purity - the true value plausibly lies below PURITY_FLOOR itself, a
    range CNVs aren't reliably callable in anyway, so searching there isn't
    meaningful. log2ratio is still evaluated, though, assuming exactly
    PURITY_FLOOR (the least purity credited), rather than ignored outright.
    At neutral CN this reduces to raw_tc unchanged by construction (nothing
    to correct for); at a real gain/loss it's a genuine correction computed
    at that assumed purity. Nothing here is artificially floored to
    PURITY_FLOOR on the way out - a neutral 4% VAF is reported as 4%, not
    pushed up to 10%, and a gain/loss correction is reported as whatever it
    actually computes to, which can land above or below PURITY_FLOOR.

    param vaf: observed variant allele fraction
    param log2ratio: local CNVkit log2 copy ratio at the variant's position,
                      or None if uncovered (falls back to the neutral model)
    param normal_cn: copy number of this locus in a normal cell (2 for
                      autosomes/female chrX, 1 for chrX in an inferred male)
    param log2_noise_floor: |log2ratio| below this is treated as neutral,
                             same as log2ratio being None - see
                             LOG2_NOISE_FLOOR for why this exists
    param purity_floor: see PURITY_FLOOR
    return: (raw_tc, adjusted_tc, cn_t, mutant_copies)
            raw_tc is always vaf * 2 (today's behaviour, kept for reference),
              never floored
            adjusted_tc is the copy-number-corrected purity estimate,
              clamped to [0, 1] - computed from the continuous solve, not
              from the rounded cn_t below. Not floored to PURITY_FLOOR: it
              can read below that (raw_tc unchanged at neutral CN below the
              floor) or above it (a real gain/loss correction evaluated at
              PURITY_FLOOR can land either side)
            cn_t is a purely informational nearest-integer label for the
              local copy number found by the solve - NOT what adjusted_tc/m
              were computed from (that's the continuous value)
            mutant_copies is the continuous m that was actually used to
              compute adjusted_tc (mutually consistent with it), which may
              not exactly match cn_t's rounded integer
    """
    raw_tc = vaf * 2

    if log2ratio is None or abs(log2ratio) < log2_noise_floor:
        cn_t, m, adjusted_tc = _cn_and_adjusted_tc(vaf, None, normal_cn, 1.0)
        return raw_tc, adjusted_tc, round(cn_t), m

    if raw_tc < purity_floor:
        # Too low to trust searching for a self-consistent purity (the true
        # value likely lies below purity_floor itself, outside the range
        # we're willing to search or report). But log2ratio is still real
        # data - rather than ignoring it outright, evaluate the correction
        # assuming exactly purity_floor (the least purity we're willing to
        # credit), and report whatever that implies, with no further floor
        # on the result: at neutral CN this naturally reduces to raw_tc
        # unchanged (nothing to correct for), and at a real gain/loss it's a
        # genuine correction computed at that assumed purity - not an
        # artificial push up to purity_floor regardless of the CN signal.
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
        # No sign change over [purity_floor, 1] - no self-consistent purity
        # exists under this model. Fall back to the uncorrected estimate.
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

    # Report adjusted_tc/m straight from the continuous solve - mutually
    # consistent and accurate. cn_t is a separate, purely informational
    # nearest-integer label; re-deriving m/adjusted_tc FROM that rounded
    # integer would reintroduce the same instability solving continuously
    # was meant to avoid (e.g. a continuous CN of 0.46 is a well-behaved
    # partial-loss state, but snaps to the degenerate cn_t=0 homozygous-
    # deletion regime if naively rounded first and recomputed from there).
    cn_t_continuous, m, adjusted_tc = _cn_and_adjusted_tc(vaf, log2ratio, normal_cn, p_star, purity_floor)
    cn_t = round(cn_t_continuous)

    return raw_tc, adjusted_tc, cn_t, m


def read_snv_vcf_and_find_max_af(input_snv_vcf, filter_dict):
    snv_vcf = pysam.VariantFile(input_snv_vcf)

    best_variant = []

    # Create VEP annotation header dict
    vep_fields = {}
    for record in snv_vcf.header.records:
        if record.type == "INFO":
            if record['ID'] == "CSQ":
                vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split('">')[0].split("|"))}

    # Iterate over the VCF file
    for record in snv_vcf.fetch():
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


def write_tc(output_tc, raw_tc, adjusted_tc):
    '''
    Write the raw and copy-number-adjusted SNV-based TC to file. Returns the
    output string to simplify unit testing.

    param output_tc: output filename
    param raw_tc: TC assuming a plain diploid heterozygous locus (VAF * 2)
    param adjusted_tc: TC corrected for local copy number
    return: output string used by the unit testing
    '''
    output = open(output_tc, "w")
    output.write("Percentage ctDNA based on SNV data (raw)\tPercentage ctDNA based on SNV data (adjusted)\n")
    output.write(f"{raw_tc*100:.1f}%\t{adjusted_tc*100:.1f}%\n")
    output.close()
    return f"{raw_tc*100:.1f}%\t{adjusted_tc*100:.1f}%\n"


# Writes additional info to file
def write_ctDNA_fraction_info(output_file_name, snv_info_list):
    '''
    Write additional info to file regarding the SNV candidates used to
    estimate TC, including the copy-number-correction details for each.

    param output_file_name: output filename
    param snv_info_list: [[raw_tc, adjusted_tc, cn_t, mutant_copies, normal_cn, VCF_record], ...]
    return: None
    '''
    output = open(output_file_name, "w")
    output.write("SNVs passing all filtering\n")
    output.write("raw_pct\tadjusted_pct\tlocal_CN_t\tassumed_mutant_copies\tnormal_CN_used\tVCF_record\n")
    for raw_tc, adjusted_tc, cn_t, m, normal_cn, record_str in snv_info_list:
        output.write(f"{raw_tc*100:.1f}%\t{adjusted_tc*100:.1f}%\t{cn_t}\t{m:.2f}\t{normal_cn}\t{record_str}")
    output.close()


if __name__ == "__main__":
    input_vcf = snakemake.input.vcf
    input_cnvkit_cns = snakemake.input.cnvkit_cns
    input_cnvkit_cnr = snakemake.input.cnvkit_cnr
    output_ctDNA_fraction = snakemake.output.ctDNA_fraction
    output_ctDNA_fraction_info = snakemake.output.ctDNA_fraction_info

    callers = snakemake.params.callers
    if isinstance(callers, list):
        callers = callers[0]

    # Building filter_dict from snakemake.params
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
        "SN": ["min", snakemake.params.min_sn],
        "AF": ["max", snakemake.params.max_af],
        "MAX_AF": ["max", snakemake.params.max_gnomad_af],
        "Consequence": ["exact", snakemake.params.excluded_consequences],
        "CHIP_genes": ["", snakemake.params.chip_genes],
        "Other": ["", []]
    }

    # Read SNVs from vcf, sorted by raw VAF descending (unchanged selection criterion)
    snv_candidates = read_snv_vcf_and_find_max_af(input_vcf, filter_dict)

    # Local copy-number lookup (CNVkit, replacing Jumble) and sex inference for chrX handling
    cns_dict = read_cnvkit_cns(input_cnvkit_cns)
    inferred_sex = infer_sex_from_cnr(input_cnvkit_cnr)

    snv_info_list = []
    for af, chrom, pos, record_str in snv_candidates:
        normal_cn = 1 if (chrom.lstrip("chr") == "X" and inferred_sex == "male") else 2
        log2ratio = lookup_local_log2(cns_dict, chrom, pos)
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(
            af, log2ratio, normal_cn, snakemake.params.log2_noise_floor, snakemake.params.purity_floor
        )
        snv_info_list.append([raw_tc, adjusted_tc, cn_t, m, normal_cn, record_str])

    if snv_info_list:
        raw_tc, adjusted_tc = snv_info_list[0][0], snv_info_list[0][1]
    else:
        raw_tc, adjusted_tc = 0, 0

    write_tc(output_ctDNA_fraction, raw_tc, adjusted_tc)
    write_ctDNA_fraction_info(output_ctDNA_fraction_info, snv_info_list)
