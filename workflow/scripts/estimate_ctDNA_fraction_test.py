import tempfile
import unittest


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.vcf = ".tests/units/estimate_ctDNA_fraction/sample1_T.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.filter.snv_hard_filter_umi.codon_snvs.sorted.vep_annotated.qci.vcf"  # noqa
        self.cnvkit_cns = ".tests/units/estimate_ctDNA_fraction/sample1_T.cns"
        self.cnvkit_cnr = ".tests/units/estimate_ctDNA_fraction/sample1_T.cnr"
        self.ctDNA_fraction = ".tests/units/estimate_ctDNA_fraction/sample1.ctDNA_fraction.tsv"

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def test_read_cnvkit_cns(self):
        from estimate_ctDNA_fraction import read_cnvkit_cns

        cns_dict = read_cnvkit_cns(self.cnvkit_cns)

        try:
            self.assertEqual([0, 30000000, 0.5849625007211562], cns_dict["chr12"][0])
            self.assertEqual([30000000, 140000000, 0.0], cns_dict["chr12"][1])
            self.assertEqual([0, 200000000, -1.0], cns_dict["chr5"][0])
        except AssertionError as e:
            print(f"Failed reading CNVkit cns. {cns_dict}")
            raise e

    def test_lookup_local_log2(self):
        from estimate_ctDNA_fraction import read_cnvkit_cns, lookup_local_log2

        cns_dict = read_cnvkit_cns(self.cnvkit_cns)

        # Position covered by a gained segment
        log2 = lookup_local_log2(cns_dict, "chr12", 25398284)
        try:
            self.assertEqual(0.5849625007211562, log2)
        except AssertionError as e:
            print(f"Failed looking up covered position. {log2}")
            raise e

        # Position not covered by any segment (chromosome not in file)
        log2 = lookup_local_log2(cns_dict, "chr21", 1000)
        try:
            self.assertIsNone(log2)
        except AssertionError as e:
            print(f"Failed looking up uncovered chromosome. {log2}")
            raise e

    def test_infer_sex_from_cnr(self):
        from estimate_ctDNA_fraction import infer_sex_from_cnr

        # sample1_T.cnr has chrX bins ~1 copy relative to autosomes -> male
        sex = infer_sex_from_cnr(self.cnvkit_cnr)
        try:
            self.assertEqual("male", sex)
        except AssertionError as e:
            print(f"Failed inferring male sex. {sex}")
            raise e

        # A sample with chrX at the same level as autosomes -> female (also the safe default)
        female_cnr = f"{self.tempdir}/female.cnr"
        with open(female_cnr, "w") as f:
            f.write("chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n")
            for i in range(12):
                f.write(f"chr1\t{i*1000}\t{i*1000+500}\t-\t0.0\t500\t1.0\n")
            for i in range(12):
                f.write(f"chrX\t{i*1000}\t{i*1000+500}\t-\t0.0\t500\t1.0\n")

        sex = infer_sex_from_cnr(female_cnr)
        try:
            self.assertEqual("female", sex)
        except AssertionError as e:
            print(f"Failed inferring female sex. {sex}")
            raise e

        # Too little chrX data to tell -> default to the safe, non-destructive "female" fallback
        sparse_cnr = f"{self.tempdir}/sparse.cnr"
        with open(sparse_cnr, "w") as f:
            f.write("chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n")
            f.write("chr1\t0\t500\t-\t0.0\t500\t1.0\n")
            f.write("chrX\t0\t500\t-\t-1.0\t500\t1.0\n")

        sex = infer_sex_from_cnr(sparse_cnr)
        try:
            self.assertEqual("female", sex)
        except AssertionError as e:
            print(f"Failed defaulting on sparse chrX data. {sex}")
            raise e

    def test_correct_vaf_for_copy_number(self):
        from estimate_ctDNA_fraction import correct_vaf_for_copy_number

        # Neutral boundary: must reduce exactly to today's VAF*2 behaviour
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.2, 0.0, 2)
        try:
            self.assertEqual(0.4, raw_tc)
            self.assertEqual(0.4, adjusted_tc)
            self.assertEqual(2, cn_t)
            self.assertEqual(1, m)
        except AssertionError as e:
            print(f"Failed neutral boundary check. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # No CNVkit coverage falls back to the neutral model
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.2, None, 2)
        try:
            self.assertEqual(0.4, raw_tc)
            self.assertEqual(0.4, adjusted_tc)
            self.assertEqual(2, cn_t)
            self.assertEqual(1, m)
        except AssertionError as e:
            print(f"Failed no-coverage fallback. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Feasible gain: vaf/log2 constructed so a true purity=0.3, CN_t=4 pair is
        # genuinely self-consistent - the continuous solve should recover both
        # almost exactly. cn_t is rounded for reporting; m is the continuous
        # value actually used to compute adjusted_tc (see correct_vaf_for_copy_number's
        # docstring for why m is deliberately not re-derived from the rounded cn_t).
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.34615384615384615, 0.3785116232537298, 2)
        try:
            self.assertAlmostEqual(0.6923076923076923, raw_tc, places=6)
            self.assertEqual(4, cn_t)
            self.assertAlmostEqual(3.0, m, places=4)
            self.assertAlmostEqual(0.3, adjusted_tc, places=4)
        except AssertionError as e:
            print(f"Failed feasible gain correction. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Feasible loss: vaf/log2 constructed so a true purity=0.3, CN_t=1 pair
        # (a simple het loss/LOH) is genuinely self-consistent.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.17647058823529413, -0.23446525363702297, 2)
        try:
            self.assertAlmostEqual(0.35294117647058826, raw_tc, places=6)
            self.assertEqual(1, cn_t)
            self.assertAlmostEqual(1.0, m, places=4)
            self.assertAlmostEqual(0.3, adjusted_tc, places=4)
        except AssertionError as e:
            print(f"Failed feasible loss correction. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Infeasible gain: log2=1.0 combined with vaf=0.3 has no self-consistent
        # purity anywhere in [PURITY_FLOOR, 1] under this model (h(p) is
        # single-signed across the whole range - previously this diverged to
        # cn_t=29 under naive fixed-point iteration). Must fall back to the
        # uncorrected neutral estimate rather than extrapolate to a boundary.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.3, 1.0, 2)
        try:
            self.assertEqual(0.6, raw_tc)
            self.assertEqual(0.6, adjusted_tc)
            self.assertEqual(2, cn_t)
            self.assertEqual(1, m)
        except AssertionError as e:
            print(f"Failed infeasible-gain fallback. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Below LOG2_NOISE_FLOOR: log2=0.02 is well within ordinary CNVkit
        # segment-level noise, not a real signal - must be treated as neutral
        # (same as log2ratio=None), regardless of how low raw_tc is. Without
        # this gate, this exact input used to manufacture a spurious
        # "correction" down to ~3.1% purely from noise crossing a rounding
        # boundary near PURITY_FLOOR.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.02, 0.02, 2)
        try:
            self.assertEqual(0.04, raw_tc)
            self.assertEqual(0.04, adjusted_tc)
            self.assertEqual(2, cn_t)
            self.assertEqual(1.0, m)
        except AssertionError as e:
            print(f"Failed log2-noise-floor gate. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Below PURITY_FLOOR but ABOVE LOG2_NOISE_FLOOR: too low to search for
        # a self-consistent purity, but log2ratio (0.15, a real signal) is
        # still evaluated - at PURITY_FLOOR itself (the least purity credited
        # down here), not ignored outright. Nothing is floored on the way
        # out: the correction pulls adjusted_tc below raw_tc, same direction
        # as it would above the floor - not pushed up to 10%.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.02, 0.15, 2)
        try:
            self.assertEqual(0.04, raw_tc)
            self.assertAlmostEqual(0.012708249882206802, adjusted_tc, places=4)
            self.assertEqual(4, cn_t)
            self.assertAlmostEqual(3.1913894413569004, m, places=4)
        except AssertionError as e:
            print(f"Failed below-purity-floor real-signal case. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # log2_noise_floor is configurable: lowering it should let a
        # previously-filtered signal (log2=0.02) through to real correction.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.02, 0.02, 2, log2_noise_floor=0.01)
        try:
            self.assertEqual(0.04, raw_tc)
            self.assertNotEqual(0.04, adjusted_tc)
        except AssertionError as e:
            print(f"Failed configurable log2_noise_floor. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Below PURITY_FLOOR with genuinely neutral CN: must reduce exactly to
        # raw_tc unchanged - no CN signal means nothing to correct for, so
        # there's no reason to report anything other than the direct VAF*2
        # value, floor or no floor.
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.04, 0.0, 2)
        try:
            self.assertEqual(0.08, raw_tc)
            self.assertAlmostEqual(0.08, adjusted_tc, places=6)
            self.assertEqual(2, cn_t)
            self.assertAlmostEqual(1.0, m, places=6)
        except AssertionError as e:
            print(f"Failed below-purity-floor neutral case. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

        # Hemizygous chrX in an inferred male (normal_cn=1): the raw VAF is directly the tumor
        # fraction, no doubling - this is the fix for the STAG2-style chrX bug
        raw_tc, adjusted_tc, cn_t, m = correct_vaf_for_copy_number(0.3, 0.0, 1)
        try:
            self.assertEqual(0.6, raw_tc)
            self.assertEqual(1, cn_t)
            self.assertAlmostEqual(1.0, m, places=6)
            self.assertAlmostEqual(0.3, adjusted_tc, places=6)
        except AssertionError as e:
            print(f"Failed hemizygous chrX correction. {raw_tc} {adjusted_tc} {cn_t} {m}")
            raise e

    def test_read_snv_vcf_and_find_max_af(self):
        from estimate_ctDNA_fraction import read_snv_vcf_and_find_max_af

        filter_dict = {
            "PositionNrSD": ["min", 20],
            "PanelMedian": ["max", 0.002],
            "Artifact": ["max", 0],
            "CALLERS": ["present", "vardict"],
            "MQ": ["min", 40],
            "MSI": ["max", 4],
            "NM": ["max", 3.0],
            "ODDRATIO": ["max", 1.5],
            "PMEAN": ["min", 25],
            "QUAL": ["min", 40],
            "SBF": ["min", 0.05],
            "SN": ["min", 100],
            "AF": ["max", 0.4],
            "MAX_AF": ["max", 0.0002],
            "Consequence": ["exact", [
                                    "synonymous_variant",
                                    "5_prime_UTR_variant",
                                    "3_prime_UTR_variant",
                                    "non_coding_transcript_exon_variant"
                                    ]],
            "CHIP_genes": ["", ["DNMT3A", "TET2", "ASXL1", "PPM1D"]],
            "Other": ["", []]
        }

        best_variant = read_snv_vcf_and_find_max_af(self.vcf, filter_dict)
        print(best_variant)

        # Only the KRAS G12V hotspot variant survives: the TET2 candidate is excluded as a CHIP
        # gene, and the others fail quality/panel-median/artifact/position filters.
        test_AF = 1.41
        test_chrom = "chr12"
        test_pos = 25398284

        try:
            self.assertEqual(1, len(best_variant))
            self.assertAlmostEqual(test_AF, best_variant[0][0] * 100, places=2)
            self.assertEqual(test_chrom, best_variant[0][1])
            self.assertEqual(test_pos, best_variant[0][2])
        except AssertionError as e:
            print(f"Failed reading vcf. {test_AF} {best_variant}")
            raise e

    def test_write_tc(self):
        from estimate_ctDNA_fraction import write_tc

        tc_string = write_tc(self.ctDNA_fraction, 0.09, 0.10)

        test_tc_string = "9.0%\t10.0%\n"

        try:
            self.assertEqual(tc_string, test_tc_string)
        except AssertionError as e:
            print(f"Failed to write output vcf. {tc_string} {test_tc_string}")
            raise e
