import logging
from pysam import AlignmentFile, VariantFile
from scipy.stats import fisher_exact


log = logging.getLogger()


def is_duplex_confirmed(read, min_per_strand_reads):
    """
    A read is duplex confirmed if fgbio's CallDuplexConsensusReads merged both original
    strands (aD/bD tags) into this consensus read. Such reads are canonicalized to a fixed
    reference-strand orientation by fgbio, so their forward/reverse alignment flag no longer
    reflects independent strand-of-origin evidence and must be excluded from a strand-bias test.
    """
    try:
        a_depth = read.get_tag("aD")
        b_depth = read.get_tag("bD")
    except KeyError:
        return False
    return a_depth >= min_per_strand_reads and b_depth >= min_per_strand_reads


def classify_read_base(read, ref_pos, ref_base, alt_base):
    for query_pos, aligned_ref_pos in read.get_aligned_pairs(matches_only=True):
        if aligned_ref_pos == ref_pos:
            base = read.query_sequence[query_pos]
            if base == alt_base:
                return "alt"
            if base == ref_base:
                return "ref"
            return None
    return None


def count_nodup_strand_alleles(bam, chrom, ref_pos, ref_base, alt_base, min_per_strand_reads):
    """
    Tally ref/alt reads by alignment strand, restricted to reads that are NOT fgbio
    duplex-confirmed (see is_duplex_confirmed) -- i.e. reads whose forward/reverse flag
    still reflects genuine, uncorrupted strand-of-origin evidence.
    """
    ref_fwd = ref_rev = alt_fwd = alt_rev = 0
    for read in bam.fetch(chrom, ref_pos, ref_pos + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if is_duplex_confirmed(read, min_per_strand_reads):
            continue
        allele = classify_read_base(read, ref_pos, ref_base, alt_base)
        if allele == "ref":
            if read.is_reverse:
                ref_rev += 1
            else:
                ref_fwd += 1
        elif allele == "alt":
            if read.is_reverse:
                alt_rev += 1
            else:
                alt_fwd += 1
    return ref_fwd, ref_rev, alt_fwd, alt_rev


def compute_nodup_annotations(ref_fwd, ref_rev, alt_fwd, alt_rev, min_total_reads):
    """
    Returns a dict with whichever of ald_nodup / (sbf_nodup, oddratio_nodup) have enough
    non-duplex-confirmed read support to be meaningful; keys for unavailable values are
    simply omitted, leaving the corresponding INFO field unset (NA) on that record.
    """
    result = {}

    if alt_fwd + alt_rev >= min_total_reads:
        result["ald_nodup"] = (alt_fwd, alt_rev)

    total = ref_fwd + ref_rev + alt_fwd + alt_rev
    if total >= min_total_reads and (alt_fwd + alt_rev) > 0 and (ref_fwd + ref_rev) > 0:
        oddsratio, pvalue = fisher_exact([[ref_fwd, ref_rev], [alt_fwd, alt_rev]])
        result["sbf_nodup"] = round(pvalue, 5)
        result["oddratio_nodup"] = round(oddsratio, 5)

    return result


def annotate(in_vcf_filename, bam_filename, out_vcf_filename, min_per_strand_reads, min_total_reads):
    in_vcf = VariantFile(in_vcf_filename)
    bam = AlignmentFile(bam_filename, "rb")

    new_header = in_vcf.header
    new_header.info.add(
        "SBF_NODUP",
        "1",
        "Float",
        "Strand bias Fisher's exact test p-value computed only from reads not confirmed as fgbio "
        "duplex consensus (aD/bD), since duplex-confirmed reads are canonicalized to a fixed strand "
        "orientation and carry no independent strand-of-origin information. Unset if too few "
        "non-duplex reads cover the position.",
    )
    new_header.info.add(
        "ODDRATIO_NODUP",
        "1",
        "Float",
        "Strand bias odds ratio, computed the same way and under the same conditions as SBF_NODUP.",
    )
    new_header.info.add(
        "ALD_NODUP",
        "2",
        "Integer",
        "Forward,reverse ALT read counts computed only from reads not confirmed as fgbio duplex "
        "consensus (aD/bD), analogous to VarDict's ALD but excluding reads whose strand orientation "
        "has been canonicalized by duplex consensus calling. Unset if too few non-duplex ALT reads "
        "cover the position.",
    )

    out_vcf = VariantFile(out_vcf_filename, "w", header=new_header)

    for record in in_vcf:
        ref_base = record.ref
        alt_base = record.alts[0] if record.alts else None
        if alt_base is not None and len(ref_base) == 1 and len(alt_base) == 1:
            ref_fwd, ref_rev, alt_fwd, alt_rev = count_nodup_strand_alleles(
                bam, record.chrom, record.pos - 1, ref_base, alt_base, min_per_strand_reads
            )
            annotations = compute_nodup_annotations(ref_fwd, ref_rev, alt_fwd, alt_rev, min_total_reads)
            if "ald_nodup" in annotations:
                record.info["ALD_NODUP"] = annotations["ald_nodup"]
            if "sbf_nodup" in annotations:
                record.info["SBF_NODUP"] = annotations["sbf_nodup"]
                record.info["ODDRATIO_NODUP"] = annotations["oddratio_nodup"]
        out_vcf.write(record)

    out_vcf.close()
    in_vcf.close()
    bam.close()


if __name__ == "__main__":
    annotate(
        snakemake.input.vcf,
        snakemake.input.bam,
        snakemake.output.vcf,
        int(snakemake.params.min_per_strand_reads),
        int(snakemake.params.min_total_reads),
    )
