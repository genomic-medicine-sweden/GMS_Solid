import os
import sys
import tempfile
import unittest

import pysam
from scipy.stats import fisher_exact

# Add the script directory to sys.path to import it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from duplex_strand_bias import (  # noqa: E402
    annotate,
    classify_read_base,
    compute_nodup_annotations,
    count_nodup_strand_alleles,
    is_duplex_confirmed,
)

HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]})
READ_START = 100
READ_LENGTH = 10
SNV_POS = 105  # 0-based; VCF position 106


def make_read(name, base_at_snv="A", reverse=False, a_depth=None, b_depth=None):
    """10 bp read covering 100-109 with the given base at SNV_POS and optional fgbio duplex tags."""
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    seq = ["A"] * READ_LENGTH
    seq[SNV_POS - READ_START] = base_at_snv
    read.query_sequence = "".join(seq)
    read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read.reference_id = 0
    read.reference_start = READ_START
    read.mapping_quality = 60
    read.cigarstring = f"{READ_LENGTH}M"
    read.is_reverse = reverse
    if a_depth is not None:
        read.set_tag("aD", a_depth)
    if b_depth is not None:
        read.set_tag("bD", b_depth)
    return read


class TestIsDuplexConfirmed(unittest.TestCase):
    def test_both_strands_supported(self):
        self.assertTrue(is_duplex_confirmed(make_read("r", a_depth=2, b_depth=1), 1))

    def test_one_strand_below_threshold(self):
        self.assertFalse(is_duplex_confirmed(make_read("r", a_depth=2, b_depth=1), 2))

    def test_missing_tags(self):
        self.assertFalse(is_duplex_confirmed(make_read("r"), 1))
        self.assertFalse(is_duplex_confirmed(make_read("r", a_depth=3), 1))


class TestClassifyReadBase(unittest.TestCase):
    def test_ref_alt_and_other(self):
        self.assertEqual(classify_read_base(make_read("r", "A"), SNV_POS, "A", "T"), "ref")
        self.assertEqual(classify_read_base(make_read("r", "T"), SNV_POS, "A", "T"), "alt")
        self.assertIsNone(classify_read_base(make_read("r", "G"), SNV_POS, "A", "T"))

    def test_position_not_covered(self):
        self.assertIsNone(classify_read_base(make_read("r"), 500, "A", "T"))


class TestComputeNodupAnnotations(unittest.TestCase):
    def test_all_annotations(self):
        result = compute_nodup_annotations(3, 2, 4, 0, min_total_reads=3)
        oddsratio, pvalue = fisher_exact([[3, 2], [4, 0]])
        self.assertEqual(result["ald_nodup"], (4, 0))
        self.assertEqual(result["sbf_nodup"], round(pvalue, 5))
        self.assertEqual(result["oddratio_nodup"], round(oddsratio, 5))

    def test_too_few_alt_reads_for_ald(self):
        result = compute_nodup_annotations(10, 10, 1, 1, min_total_reads=3)
        self.assertNotIn("ald_nodup", result)
        self.assertIn("sbf_nodup", result)

    def test_too_few_reads_overall(self):
        self.assertEqual(compute_nodup_annotations(1, 0, 1, 0, min_total_reads=3), {})

    def test_no_ref_reads_skips_fisher(self):
        result = compute_nodup_annotations(0, 0, 3, 2, min_total_reads=3)
        self.assertEqual(result, {"ald_nodup": (3, 2)})


class TestAnnotate(unittest.TestCase):
    """End-to-end: a small indexed BAM and VCF, with duplex reads that must not count."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.bam = os.path.join(self.temp_dir.name, "sample.bam")
        self.in_vcf = os.path.join(self.temp_dir.name, "in.vcf")
        self.out_vcf = os.path.join(self.temp_dir.name, "out.vcf")

        reads = (
            [make_read(f"ref_fwd{i}", "A") for i in range(3)]
            + [make_read(f"ref_rev{i}", "A", reverse=True) for i in range(2)]
            + [make_read(f"alt_fwd{i}", "T") for i in range(4)]
            # Duplex-confirmed ALT reads, all reverse: would create strand bias if counted
            + [make_read(f"alt_duplex{i}", "T", reverse=True, a_depth=2, b_depth=2) for i in range(5)]
        )
        with pysam.AlignmentFile(self.bam, "wb", header=HEADER) as out:
            for read in reads:
                out.write(read)
        pysam.index(self.bam)

        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add("chr1", length=1000)
        with pysam.VariantFile(self.in_vcf, "w", header=vcf_header) as out:
            for pos, ref, alt in [(SNV_POS + 1, "A", "T"), (SNV_POS + 1, "AA", "A"), (500, "A", "T")]:
                out.write(out.new_record(contig="chr1", start=pos - 1, alleles=(ref, alt)))

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_count_excludes_duplex_reads(self):
        with pysam.AlignmentFile(self.bam) as bam:
            counts = count_nodup_strand_alleles(bam, "chr1", SNV_POS, "A", "T", min_per_strand_reads=1)
        self.assertEqual(counts, (3, 2, 4, 0))

    def test_annotate(self):
        annotate(self.in_vcf, self.bam, self.out_vcf, min_per_strand_reads=1, min_total_reads=3)
        with pysam.VariantFile(self.out_vcf) as vcf:
            for key in ("SBF_NODUP", "ODDRATIO_NODUP", "ALD_NODUP"):
                self.assertIn(key, vcf.header.info)
            snv, indel, uncovered = list(vcf)

        oddsratio, pvalue = fisher_exact([[3, 2], [4, 0]])
        self.assertEqual(snv.info["ALD_NODUP"], (4, 0))
        self.assertAlmostEqual(snv.info["SBF_NODUP"], round(pvalue, 5), places=5)
        self.assertAlmostEqual(snv.info["ODDRATIO_NODUP"], round(oddsratio, 5), places=5)

        # Indels and positions without read support are passed through without annotations
        for record in (indel, uncovered):
            for key in ("SBF_NODUP", "ODDRATIO_NODUP", "ALD_NODUP"):
                self.assertNotIn(key, record.info)

    def test_duplex_reads_counted_when_threshold_not_met(self):
        # With min_per_strand_reads above the aD/bD depths, the duplex reads are not
        # considered duplex-confirmed and their reverse ALT support is counted.
        with pysam.AlignmentFile(self.bam) as bam:
            counts = count_nodup_strand_alleles(bam, "chr1", SNV_POS, "A", "T", min_per_strand_reads=3)
        self.assertEqual(counts, (3, 2, 4, 5))


if __name__ == "__main__":
    unittest.main()
