
__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"

import pysam
import statistics
import numpy as np
import scipy.stats as stats


def read_segments(input_segments):
    """
    Read segments from CNV segmentation file and store in dict. Skips X and Y chromosome.

    param input_segments: File handle to CNV segments
    return segment_dict: a dict (with chromosomes as keys) each with a list of segments.
                Each segment has: [start position, end position, copy number, [empty list where germline AFs will be stored]]
    """
    segments = pysam.VariantFile(input_segments)

    segment_dict = {}
    for segment in segments.fetch():
        chrom = segment.chrom
        start_pos = segment.pos
        end_pos = segment.stop
        CN = segment.info["CORR_CN"]
        if chrom.endswith("X") or chrom.endswith("Y"):
            continue
        if chrom not in segment_dict:
            segment_dict[chrom] = []
        segment_dict[chrom].append([start_pos, end_pos, CN, []])
    return segment_dict


def read_germline_vcf(input_germline_vcf, segment_dict, min_germline_af):
    """
    Read germline SNPs from vcf file. Skips low (min_germline_af) and high AF (1 - min_germline_af) SNPs.
    Store SNPs in dict and also update the segment dict with germline AFs within each segment.

    param input_germline_vcf: File handle to vcf file
    param segment_dict: Dict created by the read_segments function
    param min_germline_af: Float with minimum AF to be counted as a germline.
    return segment_dict: a dict (with chromosomes as keys) each with a list of segments.
                Each segment has: [start position, end position, copy number, [list of germline AFs]]
    return germline_dict: a dict (with chromosomes as keys) each with a list of SNP-data.
                The SNP-data contains: [position, AF]
    """
    germline_vcf = pysam.VariantFile(input_germline_vcf)
    germline_dict = {}

    for record in germline_vcf.fetch():
        chrom = record.chrom
        pos = record.pos
        AF = record.info["AF"][0]
        if AF < min_germline_af or AF > (1 - min_germline_af):
            continue
        if chrom in segment_dict:
            i = 0
            for segment in segment_dict[chrom]:
                if pos >= segment[0] and pos <= segment[1]:
                    segment_dict[chrom][i][3].append(AF)
                i += 1
        if chrom not in germline_dict:
            germline_dict[chrom] = []
        germline_dict[chrom].append([pos, AF])
    return segment_dict, germline_dict


def read_bedfile(filter_regions_dict, in_bed_filename):
    """
    Read a bed file and add these regions to the filter_regions_dict.

    param filter_regions_dict: A dict (with chromosomes as keys) each with a [list of regions: [start position, end position]]
    param in_bed_filename: File name of the file
    return filter_regions_dict: Updated dict with additional regions
    """
    in_bed = open(in_bed_filename)
    for line in in_bed:
        columns = line.strip().split("\t")
        chrom = columns[0]
        if chrom not in filter_regions_dict:
            filter_regions_dict[chrom] = []
        filter_regions_dict[chrom].append([int(columns[1]), int(columns[2])])
    in_bed.close()
    return filter_regions_dict


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
                        ("Hotspot" in record.info and (
                            record.info["Hotspot"] == "3-check" or
                            record.info["Hotspot"] == "1-hotspot"
                        ))
                        ):
                    filtered = True
            elif filter == "CHIP_genes":
                if vep_dict["SYMBOL"] in filter_dict[filter][1]:
                    filtered = True

        if not filtered:
            best_variant.append([record.info["AF"][0], str(record)])

    if not best_variant:
        return 0, []

    best_variant.sort(key=lambda x: x[0], reverse=True)
    return best_variant[0][0] * 2, best_variant


def baf_to_tc(abs_value_seg_median, CN, CN_list, median_noise_level):
    '''
    Translates the baf signal size of a region into a TC value.
    The baf signal is reduced by the background noise level.
    Then, depending on the nature of the CNA, the baf signal is transformed into a TC value.

    param abs_value_seg_median: Float with the median absolute baf signal (signal = difference from 50% AF)
    param CN: Copy number of the signal
    param CN_list: List of copy numbers of all segments with baf signal
    param median_noise_level: Float with the median baf noise level based on germline AF in copy neutral segments
    return: [TC, CNA-type]
    '''
    # Remove the medium noise level
    abs_value_seg_median -= median_noise_level
    # Get highest and lowest difference to normal copy number for a segment
    min_CN_diff = 0.00001
    max_CN_diff = 0.00001
    if len(CN_list) > 0:
        min_CN_diff = max(0.00001, 2.0 - min(CN_list))
        max_CN_diff = max(0.00001, max(CN_list) - 2.0)
    # If clear aneuploidy (only based on copy number)
    if len(CN_list) > 0 and (min(CN_list) < 1.85 or max(CN_list) > 2.15):
        # Deletion if copy number of at least 40% of the lowest copy number segment
        if (2.0 - CN) / min_CN_diff > 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Del"]
        # Duplication (not reported) if copy number of at least 40% of the highest copy number segment
        elif (CN - 2.0) / max_CN_diff > 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.333), "Dup"]
        # Copy neutral LoH (Reported if no deletion) if copy number looks normal
        elif (2.0 - CN) / min_CN_diff < 0.4 and (CN - 2.0) / max_CN_diff < 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.667), "CNLoH"]
        # Unknown, assuming deletion (not reported)
        else:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Unknown"]
    # If no clear CNA signal, assume deletion for TC adjustment model for all segments
    else:
        return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Del"]


def test_if_signal_in_segment(segment, abs_value_seg_median, median_noise_level, vaf_baseline):
    '''
    Test if the baf signal in a segment is due to a CNA or just noise.
    Uses density curves of baf-values to determine if there are two maxima centered with a minima around the VAF baseline (0.5)

    param segment: [list of germline AF in segment]
    param abs_value_seg_median: Float with the median absolute baf signal (signal = difference from 50% AF)
    param median_noise_level: Float with the median baf noise level based on germline AF in copy neutral segments
    param vaf_baseline: Float that describes where the baf baseline are.
                        In theory it should be 0.5 but in practice it is usually slightly lower (0.48).

    return True (the segment has CNA) if all of these are true (otherwise False):
        * the signal is higher than the background noise
        * the minimum VAF is not the same as any maxima (only one peak)
        * the minimum VAF is between 0.4 and 0.6 (peaks should be centered around the VAF baseline)
        * the maxima ratio (max1 / max2) is between 0.5 and 2 (there should be similar number of
                                                               germline SNPs on both sides of the VAF baseline)
        * both the minimum maximum ratio (min / max1 and min / max2) is lower than 0.85
    '''
    # Get VAF for max density on the lower AF half of the BAF plot
    kde1 = stats.gaussian_kde(segment)
    kde1.set_bandwidth(bw_method=kde1.factor / 2.)
    x = np.linspace(0, vaf_baseline, 200)
    kde1 = kde1(x)
    max1 = x[kde1.argmax()]

    # Get VAF for max density on the upper AF half of the BAF plot
    kde2 = stats.gaussian_kde(segment)
    kde2.set_bandwidth(bw_method=kde2.factor / 2.)
    x = np.linspace(vaf_baseline, 1, 200)
    kde2 = kde2(x)
    max2 = x[kde2.argmax()]

    # Get VAF for min density between the two maxima
    kde3 = stats.gaussian_kde(segment)
    kde3.set_bandwidth(bw_method=kde3.factor / 2.)
    x = np.linspace(max1, max2, 200)
    kde3 = kde3(x)
    minimum = x[kde3.argmin()]

    # Get the density for the maxima and minimum
    density_max1 = kde1[kde1.argmax()]
    density_max2 = kde2[kde2.argmax()]
    density_min = kde3[kde3.argmin()]

    if abs_value_seg_median > median_noise_level:
        if minimum != max1 and minimum != max2 and minimum > 0.4 and minimum < 0.6:
            if (
                density_max1 / density_max2 < 2 and
                density_max1 / density_max2 > 0.5 and
                density_min / density_max1 < 0.85 and
                density_min / density_max2 < 0.85
            ):
                return True
    return False


def calculate_cnv_tc(segment_dict_AF, min_nr_SNPs_per_segment, vaf_baseline, min_segment_length):
    """
    First calculate BAF noise levels for segments without CNA (median_noise_level)
    Then calculate TC for all segments that have signal (based on function test_if_signal_in_segment)
    Return the highest TC based segments with deletions and LoH.
    Also return a list of segments with signal and their estimated TC (seg_list)

    param segment_dict_AF: Dict created by the read_segments function and updated by the read_germline_vcf function
    param min_nr_SNPs_per_segment: Number of germline SNPs in segment needed to do a reliable baf signal test
    param vaf_baseline: Float that describes where the baf baseline are.
                        In theory it should be 0.5 but in practice it is usually slightly lower (0.48).
    min_segment_length: Minimal length of the segment (in bases) (avoid small noisy segments)
    return max_tc: The maximum TC
    return seg_list: [[TC, [segment], chrom], CNA-type]
                     segment = [start position, end position, copy number, [list of germline AFs]]
    """
    CN_signal_list = []
    noise_level = []
    for chrom in segment_dict_AF:
        for segment in segment_dict_AF[chrom]:
            if len(segment[3]) > min_nr_SNPs_per_segment:
                # Check if CNV segment has VAF signal (using 0 noise levels) and save in CN_signal_list for later reference.
                # If not significant save all background VAF signals in noise_level.
                if not test_if_signal_in_segment(segment[3], 1, 0, vaf_baseline):
                    for AF in segment[3]:
                        noise_level.append(abs(AF-vaf_baseline))
                else:
                    CN_signal_list.append(segment[2])
    # Calculated the median noise level
    if len(noise_level) > 0:
        median_noise_level = statistics.median(noise_level)
    else:
        median_noise_level = 0

    # Iterate through segments to find final CNAs
    tc_dict = {"Del": [], "Dup": [], "CNLoH": [], "Unknown": []}
    for chrom in segment_dict_AF:
        for segment in segment_dict_AF[chrom]:
            data_points_seg = len(segment[3])
            segment_length = segment[1] - segment[0]
            # Only use segments that have sufficient number of SNPs and segment length (in bases)
            if data_points_seg > min_nr_SNPs_per_segment and segment_length > min_segment_length:
                i = 0
                short = True
                # For extra long segments, split them up into 100 datapoint segments to find smaller CNAs
                while i + min_nr_SNPs_per_segment < data_points_seg:
                    seg = segment[3]
                    if data_points_seg > 100:
                        seg = segment[3][i:i+100]
                        short = False
                    abs_value_seg = []
                    for AF in seg:
                        abs_value_seg.append(abs(AF-vaf_baseline))
                    abs_value_seg_median = statistics.median(abs_value_seg)
                    # Check if signal is found and is higher than the background noise level
                    if test_if_signal_in_segment(seg, abs_value_seg_median, median_noise_level, vaf_baseline):
                        # Calculate TC based on the median separation in BAF around the BAF-baseline (~50%)
                        tc_seg = baf_to_tc(abs_value_seg_median, segment[2], CN_signal_list, median_noise_level)
                        tc_dict[tc_seg[1]].append([tc_seg[0], segment, chrom])
                    i += 50
                    if short:
                        break
    # Report highest TC based on CNAs with deletions.
    # If none are found report highest copy neutral LoH CNA.
    max_tc = 0
    found_del = False
    seg_list = []
    for seg_info in tc_dict["Del"]:
        seg_list.append([seg_info, "Deletion"])
        if seg_info[0] > max_tc:
            max_tc = seg_info[0]
            found_del = True
    for seg_info in tc_dict["CNLoH"]:
        seg_list.append([seg_info, "CNLoH"])
        if found_del and seg_info[0] > max_tc:
            max_tc = seg_info[0]
    return max_tc, seg_list


def write_tc(output_tc, tc_cnv, tc_snv):
    '''
    Write CNV and SNV based TC to file. Return output string to simplify unit testing.

    param output_tc: output filename
    param tc_cnv: TC based on CNAs
    param tc_snv: TC based on somatic SNVs
    return: output string used by the unit testing
    '''
    output = open(output_tc, "w")
    output.write("Percentage ctDNA based on CNV data\tPercentage ctDNA based on SNV data\n")
    output.write(f"{tc_cnv*100:.1f}%\t{tc_snv*100:.1f}%\n")
    output.close()
    return f"{tc_cnv*100:.1f}%\t{tc_snv*100:.1f}%\n"


# Writes additional info to file
def write_ctDNA_fraction_info(output_file_name, seg_list, snv_list):
    '''
    Write additional info to file regarding the CNA and SNV candidates used to estimate TC.

    param output_file_name: output filename
    param seg_list: [list of segments: [ctDNA_percentage, CNV_type, chromosome, start position, end position]]
    param snv_list: [list of SNVs (the entire row in the vcf)]
    return: output string used by the unit testing
    '''
    output = open(output_file_name, "w")
    output.write("ctDNA_percentage\tCNV_type\tchromosome\tstart_pos\tend_pos\n")
    for seg in seg_list:
        output.write(f"{seg[0][0]*100:.1f}%\t{seg[1]}\t{seg[0][2]}\t{seg[0][1][0]}\t{seg[0][1][1]}\n")
    output.write("\nSNVs passing all filtering\n")
    output.write("ctDNA_percentage\tVCF_record\n")
    for snv in snv_list:
        output.write(f"{snv[0]*2*100:.1f}%\t{snv[1]}")
    output.close()


if __name__ == "__main__":
    input_segments = snakemake.input.segments
    input_germline_vcf = snakemake.input.germline_vcf
    input_vcf = snakemake.input.vcf
    output_ctDNA_fraction = snakemake.output.ctDNA_fraction
    output_ctDNA_fraction_info = snakemake.output.ctDNA_fraction_info

    min_germline_af = float(snakemake.params.min_germline_af)
    min_nr_SNPs_per_segment = int(snakemake.params.min_nr_SNPs_per_segment)
    min_segment_length = int(snakemake.params.min_segment_length)
    vaf_baseline = float(snakemake.params.vaf_baseline)

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

    # Read CNV segments
    segment_dict = read_segments(input_segments)
    # Read germline SNPs
    segment_dict_AF, germline_dict = read_germline_vcf(input_germline_vcf, segment_dict, min_germline_af)
    # Calculate TC based on BAF germline values
    tc_cnv, seg_list = calculate_cnv_tc(segment_dict_AF, min_nr_SNPs_per_segment, vaf_baseline, min_segment_length)
    # Read SNVs from vcf and then report TC based on max VAF of somatic SNV.
    tc_snv, snv_list = read_snv_vcf_and_find_max_af(input_vcf, filter_dict)
    write_tc(output_ctDNA_fraction, tc_cnv, tc_snv)
    # Write additional info regarding which chromosomes have deletions
    # Write additional info regarding additional SNVs found
    write_ctDNA_fraction_info(output_ctDNA_fraction_info, seg_list, snv_list)
