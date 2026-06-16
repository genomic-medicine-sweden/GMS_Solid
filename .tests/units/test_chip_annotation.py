import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import pysam

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from chip_annotation import (  # noqa: E402
    get_vep_field_index,
    lookup_tabix_vcf,
    chip_gene_flag,
    chip_vaf_flag,
    chip_frag_and_sbs_flags,
    chip_hotspot_flag,
    chip_cosmic_hemato_flag,
    chip_consequence_flag,
)

# VCF header with CSQ format: Allele|Consequence|SYMBOL  (indices 0, 1, 2)
_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP. Format: Allele|Consequence|SYMBOL">\n'
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n'
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##contig=<ID=chr2,length=243199373>\n"
    "##contig=<ID=chr17,length=81195210>\n"
    "##contig=<ID=chr7,length=159138663>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)

SYMBOL_IDX = 2
CONSEQUENCE_IDX = 1
MAJOR_GENES = {"DNMT3A", "TET2", "ASXL1"}
MINOR_GENES = {"TP53", "JAK2", "SF3B1"}
LOF_GENES = {"DNMT3A", "TET2", "ASXL1"}
MISSENSE_GENES = {"TP53", "JAK2"}


def _make_vcf(rows):
    """Write a temp VCF with the standard header plus the given data rows; return the path."""
    f = tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False)
    f.write(_VCF_HEADER)
    for row in rows:
        f.write(row + "\n")
    f.close()
    return f.name


class TestGetVepFieldIndex(unittest.TestCase):
    def _header(self, csq_format):
        record = MagicMock()
        record.key = "INFO"
        record.get = lambda key, default=None: {"ID": "CSQ", "Description": f'VEP. Format: {csq_format}'}.get(key, default)
        header = MagicMock()
        header.records = [record]
        return header

    def test_symbol_found(self):
        self.assertEqual(get_vep_field_index(self._header("Allele|Consequence|SYMBOL"), "SYMBOL"), 2)

    def test_consequence_found(self):
        self.assertEqual(get_vep_field_index(self._header("Allele|Consequence|SYMBOL"), "Consequence"), 1)

    def test_field_absent(self):
        self.assertIsNone(get_vep_field_index(self._header("Allele|Consequence|SYMBOL"), "HGVSp"))

    def test_no_csq_info(self):
        record = MagicMock()
        record.key = "FILTER"  # not INFO
        header = MagicMock()
        header.records = [record]
        self.assertIsNone(get_vep_field_index(header, "SYMBOL"))


class TestLookupTabixVcf(unittest.TestCase):
    def _tabix(self, rows):
        mock = MagicMock()
        mock.fetch.return_value = iter(rows)
        return mock

    def test_exact_match(self):
        tabix = self._tabix(["chr2\t100\t.\tC\tT\t.\t.\tCNT=15"])
        self.assertEqual(lookup_tabix_vcf(tabix, "chr2", 100, "C", "T"), "CNT=15")

    def test_wrong_alt(self):
        tabix = self._tabix(["chr2\t100\t.\tC\tA\t.\t.\tCNT=15"])
        self.assertIsNone(lookup_tabix_vcf(tabix, "chr2", 100, "C", "T"))

    def test_wrong_ref(self):
        tabix = self._tabix(["chr2\t100\t.\tG\tT\t.\t.\tCNT=15"])
        self.assertIsNone(lookup_tabix_vcf(tabix, "chr2", 100, "C", "T"))

    def test_no_rows(self):
        tabix = self._tabix([])
        self.assertIsNone(lookup_tabix_vcf(tabix, "chr2", 100, "C", "T"))

    def test_fetch_raises_value_error(self):
        tabix = MagicMock()
        tabix.fetch.side_effect = ValueError("region error")
        self.assertIsNone(lookup_tabix_vcf(tabix, "chrUn", 1, "A", "T"))

    def test_no_info_column(self):
        # Row with fewer than 8 columns returns "."
        tabix = self._tabix(["chr2\t100\t.\tC\tT\t.\t."])
        self.assertEqual(lookup_tabix_vcf(tabix, "chr2", 100, "C", "T"), ".")


class TestChipGeneFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf([
            "chr2\t25234373\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|DNMT3A\tGT\t0/1",    # major
            "chr17\t7674220\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|TP53\tGT\t0/1",       # minor
            "chr7\t55242467\t.\tG\tA\t.\t.\tCSQ=A|missense_variant|EGFR\tGT\t0/1",       # non-CHIP
            "chr2\t25234374\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|TP53,G|missense_variant|DNMT3A\tGT\t0/1",  # major+minor
            "chr2\t25234375\t.\tA\tG\t.\t.\t.\tGT\t0/1",                                  # no CSQ
        ])
        cls._vcf = pysam.VariantFile(path)
        cls._recs = list(cls._vcf)
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def test_major_gene(self):
        self.assertEqual(chip_gene_flag(self._recs[0], SYMBOL_IDX, MAJOR_GENES, MINOR_GENES), 2)

    def test_minor_gene(self):
        self.assertEqual(chip_gene_flag(self._recs[1], SYMBOL_IDX, MAJOR_GENES, MINOR_GENES), 1)

    def test_non_chip_gene(self):
        self.assertEqual(chip_gene_flag(self._recs[2], SYMBOL_IDX, MAJOR_GENES, MINOR_GENES), 0)

    def test_major_beats_minor(self):
        self.assertEqual(chip_gene_flag(self._recs[3], SYMBOL_IDX, MAJOR_GENES, MINOR_GENES), 2)

    def test_no_csq(self):
        self.assertEqual(chip_gene_flag(self._recs[4], SYMBOL_IDX, MAJOR_GENES, MINOR_GENES), 0)

    def test_symbol_idx_none(self):
        self.assertEqual(chip_gene_flag(self._recs[0], None, MAJOR_GENES, MINOR_GENES), 0)


class TestChipVafFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf([
            "chr2\t1\t.\tC\tT\t.\t.\tAF=0.02\tGT\t0/1",   # 2% → tier 2
            "chr2\t2\t.\tC\tT\t.\t.\tAF=0.05\tGT\t0/1",   # 5% → tier 1
            "chr2\t3\t.\tC\tT\t.\t.\tAF=0.20\tGT\t0/1",   # 20% → tier 0
            "chr2\t4\t.\tC\tT\t.\t.\tAF=0.006\tGT\t0/1",  # 0.6% → tier 1 (near lower boundary)
            "chr2\t5\t.\tC\tT\t.\t.\t.\tGT\t0/1",         # no AF in INFO → 0
        ])
        cls._vcf = pysam.VariantFile(path)
        cls._recs = list(cls._vcf)
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def test_tier2_range(self):
        self.assertEqual(chip_vaf_flag(self._recs[0]), 2)

    def test_tier1_range(self):
        self.assertEqual(chip_vaf_flag(self._recs[1]), 1)

    def test_too_high(self):
        self.assertEqual(chip_vaf_flag(self._recs[2]), 0)

    def test_lower_boundary_tier1(self):
        self.assertEqual(chip_vaf_flag(self._recs[3]), 1)

    def test_no_af_field(self):
        self.assertEqual(chip_vaf_flag(self._recs[4]), 0)


class TestChipFragFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf([
            "chr2\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1",   # SNV record
            "chr2\t101\t.\tC\t.\t.\t.\t.\tGT\t0/1",   # no alts
        ])
        cls._vcf = pysam.VariantFile(path)
        cls._recs = list(cls._vcf)
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def _call(self, alt_lens, ref_lens, min_alt=3, abs_thr=150, ratio_thr=1.1, short_thr=145):
        with patch("chip_annotation.pileup_at_variant", return_value=(alt_lens, ref_lens, 0, 0)):
            return chip_frag_and_sbs_flags(self._recs[0], MagicMock(), min_alt, abs_thr, ratio_thr, short_thr)

    def test_no_alts_returns_none(self):
        result = chip_frag_and_sbs_flags(self._recs[1], MagicMock(), 3, 150, 1.1, 145)
        self.assertEqual(result, (None, None, None))

    def test_too_few_alt_reads(self):
        result = self._call([160], [120, 125], min_alt=5)
        self.assertEqual(result, (None, None, None))

    def test_both_criteria_met(self):
        # alt_mean=165, ref_mean=140 → abs > 150, ratio=1.18 > 1.1
        chip_frag, frag_alt, frag_ref = self._call([160, 165, 170], [135, 140, 145])
        self.assertEqual(chip_frag, 1)
        self.assertAlmostEqual(frag_alt, 165.0)
        self.assertAlmostEqual(frag_ref, 140.0)

    def test_abs_met_ratio_not(self):
        # Germline-like: alt_mean=183, ref_mean=183 → abs > 150 but ratio=1.0 < 1.1
        chip_frag, frag_alt, frag_ref = self._call([180, 183, 186], [180, 183, 186])
        self.assertEqual(chip_frag, 0)

    def test_ratio_met_abs_not(self):
        # alt_mean=147, ref_mean=127 → ratio=1.16 > 1.1 but abs < 150; above short threshold (145)
        chip_frag, frag_alt, frag_ref = self._call([144, 147, 150], [122, 127, 132])
        self.assertEqual(chip_frag, 0)

    def test_neither_criterion_met(self):
        # alt_mean=148, ref_mean=145 → abs < 150, ratio < 1.1, not below short threshold
        chip_frag, frag_alt, frag_ref = self._call([145, 148, 151], [142, 145, 148])
        self.assertEqual(chip_frag, 0)

    def test_short_fragment(self):
        # alt_mean=130 < 145 → tumour-like, chip_frag=-1
        chip_frag, frag_alt, frag_ref = self._call([125, 130, 135], [160, 165, 170])
        self.assertEqual(chip_frag, -1)

    def test_short_takes_priority_over_ratio(self):
        # alt_mean=130 < 145, even if ratio would be met → -1 wins
        chip_frag, frag_alt, frag_ref = self._call([125, 130, 135], [100, 105, 110])
        self.assertEqual(chip_frag, -1)

    def test_no_ref_reads(self):
        # No ref reads → ratio cannot be assessed → chip_frag=0
        chip_frag, frag_alt, frag_ref = self._call([160, 165, 170], [])
        self.assertIsNone(frag_ref)
        self.assertEqual(chip_frag, 0)


class TestChipHotspotFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf(["chr2\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1"])
        cls._vcf = pysam.VariantFile(path)
        cls._rec = list(cls._vcf)[0]
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def test_found_in_hotspot(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr2\t100\t.\tC\tT\t.\t.\tCNT=50"])
        self.assertEqual(chip_hotspot_flag(self._rec, tabix), 1)

    def test_not_in_hotspot(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter([])
        self.assertEqual(chip_hotspot_flag(self._rec, tabix), 0)

    def test_tabix_none(self):
        self.assertEqual(chip_hotspot_flag(self._rec, None), 0)


class TestChipCosmicHematoFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf([
            "chr2\t100\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|DNMT3A\tGT\t0/1",
            "chr17\t100\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|TP53\tGT\t0/1",
        ])
        cls._vcf = pysam.VariantFile(path)
        recs = list(cls._vcf)
        cls._rec_dnmt3a = recs[0]
        cls._rec_tp53 = recs[1]
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def test_cnt_returned(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr2\t100\t.\tC\tT\t.\t.\tCNT=42"])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, tabix, SYMBOL_IDX, set(), 20), 42)

    def test_cnt_below_min_returns_0(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr2\t100\t.\tC\tT\t.\t.\tCNT=12"])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, tabix, SYMBOL_IDX, set(), 20), 0)

    def test_cnt_at_min_returned(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr2\t100\t.\tC\tT\t.\t.\tCNT=20"])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, tabix, SYMBOL_IDX, set(), 20), 20)

    def test_no_cnt_field_returns_0(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr2\t100\t.\tC\tT\t.\t.\tGENE=DNMT3A"])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, tabix, SYMBOL_IDX, set(), 20), 0)

    def test_not_found_returns_0(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter([])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, tabix, SYMBOL_IDX, set(), 20), 0)

    def test_cosmic_tabix_none(self):
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_dnmt3a, None, SYMBOL_IDX, set(), 20), 0)

    def test_excluded_gene_returns_0(self):
        tabix = MagicMock()
        tabix.fetch.return_value = iter(["chr17\t100\t.\tC\tT\t.\t.\tCNT=456"])
        self.assertEqual(chip_cosmic_hemato_flag(self._rec_tp53, tabix, SYMBOL_IDX, {"TP53"}, 20), 0)

    def test_excluded_gene_not_queried(self):
        tabix = MagicMock()
        chip_cosmic_hemato_flag(self._rec_tp53, tabix, SYMBOL_IDX, {"TP53"}, 20)
        tabix.fetch.assert_not_called()


class TestChipConsequenceFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = _make_vcf([
            "chr2\t1\t.\tC\tT\t.\t.\tCSQ=T|frameshift_variant|DNMT3A\tGT\t0/1",    # LoF in lof_gene
            "chr17\t1\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|TP53\tGT\t0/1",        # missense in missense_gene
            "chr7\t1\t.\tC\tT\t.\t.\tCSQ=T|frameshift_variant|TP53\tGT\t0/1",       # LoF but TP53 not in lof_genes
            "chr2\t2\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|DNMT3A\tGT\t0/1",       # missense but DNMT3A not missense_gene
            "chr2\t3\t.\tC\tT\t.\t.\tCSQ=T|synonymous_variant|DNMT3A\tGT\t0/1",     # synonymous → 0
            "chr2\t4\t.\tC\tT\t.\t.\t.\tGT\t0/1",                                   # no CSQ → 0
            # compound: LoF&missense in lof_gene → 1
            "chr2\t5\t.\tC\tT\t.\t.\tCSQ=T|stop_gained&missense_variant|DNMT3A\tGT\t0/1",
        ])
        cls._vcf = pysam.VariantFile(path)
        cls._recs = list(cls._vcf)
        os.unlink(path)

    @classmethod
    def tearDownClass(cls):
        cls._vcf.close()

    def test_lof_in_lof_gene(self):
        self.assertEqual(chip_consequence_flag(self._recs[0], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 1)

    def test_missense_in_missense_gene(self):
        self.assertEqual(chip_consequence_flag(self._recs[1], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 1)

    def test_lof_not_in_lof_gene(self):
        self.assertEqual(chip_consequence_flag(self._recs[2], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 0)

    def test_missense_not_in_missense_gene(self):
        self.assertEqual(chip_consequence_flag(self._recs[3], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 0)

    def test_synonymous(self):
        self.assertEqual(chip_consequence_flag(self._recs[4], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 0)

    def test_no_csq(self):
        self.assertEqual(chip_consequence_flag(self._recs[5], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 0)

    def test_compound_consequence(self):
        self.assertEqual(chip_consequence_flag(self._recs[6], SYMBOL_IDX, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 1)

    def test_symbol_idx_none(self):
        self.assertEqual(chip_consequence_flag(self._recs[0], None, CONSEQUENCE_IDX, LOF_GENES, MISSENSE_GENES), 0)

    def test_consequence_idx_none(self):
        self.assertEqual(chip_consequence_flag(self._recs[0], SYMBOL_IDX, None, LOF_GENES, MISSENSE_GENES), 0)


if __name__ == "__main__":
    unittest.main()
