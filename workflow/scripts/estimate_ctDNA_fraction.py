
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


def correct_vaf_for_copy_number(vaf, log2ratio, normal_cn=2):
    """
    Convert a somatic VAF into a tumor-fraction estimate, correcting for
    local copy number instead of always assuming a plain diploid
    heterozygous locus (today's VAF*2 behaviour).

    Model: the mutant allele is assumed to be the one preferentially
    amplified on a gain (m = CN_t - (normal_cn - 1), i.e. exactly
    normal_cn - 1 copies stay wild-type and everything else is mutant), and
    the retained copy/copies are assumed mutant on a loss/LOH (m = CN_t).
    At CN_t == normal_cn this reduces to the plain heterozygous model (m=1),
    matching the previous VAF*2 behaviour exactly.

    param vaf: observed variant allele fraction
    param log2ratio: local CNVkit log2 copy ratio at the variant's position,
                      or None if uncovered (falls back to the neutral model)
    param normal_cn: copy number of this locus in a normal cell (2 for
                      autosomes/female chrX, 1 for chrX in an inferred male)
    return: (raw_tc, adjusted_tc, cn_t, mutant_copies)
            raw_tc is always vaf * 2 (today's behaviour, kept for reference)
            adjusted_tc is the copy-number-corrected estimate, clamped [0, 1]
            cn_t is the local copy number used (rounded to the nearest
              integer for the allele-count model)
            mutant_copies is the assumed number of mutant copies (m)
    """
    raw_tc = vaf * 2

    if log2ratio is None:
        cn_t = normal_cn
    else:
        cn_t = max(round(normal_cn * (2 ** log2ratio)), 0)

    if cn_t >= normal_cn:
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
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(af, log2ratio, normal_cn)
        snv_info_list.append([raw_tc, adjusted_tc, cn_t, m, normal_cn, record_str])

    if snv_info_list:
        raw_tc, adjusted_tc = snv_info_list[0][0], snv_info_list[0][1]
    else:
        raw_tc, adjusted_tc = 0, 0

    write_tc(output_ctDNA_fraction, raw_tc, adjusted_tc)
    write_ctDNA_fraction_info(output_ctDNA_fraction_info, snv_info_list)
