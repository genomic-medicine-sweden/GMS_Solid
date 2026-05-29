import tempfile
import os
import unittest

test_segment_dict = {
    'chr2': [[138974, 242801052, 1.9299999475479126, []]],
    'chr3': [[94842, 196812807, 1.9900000095367432, []]],
    'chr4': [[29134, 191153264, 1.9600000381469727, []]],
    'chr5': [[218403, 180914250, 2.009999990463257, []]],
}

test_segment_dict2 = {
    'chr2': [[138974, 242801052, 1.9299999475479126, []]],
    'chr3': [[94842, 196812807, 1.9900000095367432, []]],
    'chr4': [[29134, 191153264, 1.9600000381469727, []]],
    'chr5': [[218403, 180914250, 2.009999990463257, []]],
}

test_updated_dict = {
    'chr2': [[138974, 242801052, 1.9299999475479126, []]],
    'chr3': [[94842, 196812807, 1.9900000095367432, []]],
    'chr4': [[29134, 191153264, 1.9600000381469727, []]],
    'chr5': [[218403, 180914250, 2.009999990463257, []]],
}

test_segment_dict3 = {
    'chr2': [[138974, 242801052, 1.9299999475479126,
              [0.4519999921321869, 0.5094000101089478, 0.4846999943256378, 0.5286999940872192, 0.4810999929904938,
               0.42480000853538513, 0.4255000054836273, 0.4207000136375427, 0.4794999957084656, 0.4399999976158142,
               0.44769999384880066, 0.5174999833106995, 0.15629999339580536, 0.44209998846054077, 0.4932999908924103,
               0.4903999865055084]]],
    'chr3': [[94842, 196812807, 1.9900000095367432,
              [0.47519999742507935, 0.42309999465942383, 0.5034999847412109, 0.49779999256134033, 0.4634000062942505,
               0.4875999987125397, 0.4848000109195709, 0.5054000020027161, 0.4823000133037567, 0.4625999927520752,
               0.5281000137329102, 0.550000011920929, 0.46869999170303345]]],
    'chr4': [[29134, 191153264, 1.9600000381469727,
              [0.4104999899864197, 0.5218999981880188, 0.49950000643730164, 0.45730000734329224, 0.4277999997138977,
               0.6115000247955322, 0.46810001134872437, 0.49160000681877136, 0.451200008392334, 0.4607999920845032,
               0.5008000135421753, 0.5677000284194946, 0.4672999978065491, 0.49050000309944153, 0.41690000891685486,
               0.47189998626708984, 0.476500004529953, 0.4625999927520752]]],
    'chr5': [[218403, 180914250, 2.009999990463257,
              [0.5001999735832214, 0.4828999936580658, 0.490200012922287, 0.4474000036716461, 0.4496000111103058,
               0.47929999232292175, 0.534500002861023, 0.5522000193595886, 0.4812999963760376, 0.4706000089645386,
               0.5550000071525574, 0.4221000075340271, 0.4345000088214874, 0.4458000063896179]]]
}


test_segment = [0.6001999735832214, 0.3828999936580658, 0.690200012922287, 0.3474000036716461, 0.3496000111103058,
                0.37929999232292175, 0.634500002861023, 0.6522000193595886, 0.3812999963760376, 0.3706000089645386,
                0.6550000071525574, 0.3221000075340271, 0.3345000088214874, 0.3458000063896179]


test_segment_dict4 = {
    'chr2': [[138974, 242801052, 1.9299999475479126,
              [0.4519999921321869, 0.5094000101089478, 0.4846999943256378, 0.5286999940872192, 0.4810999929904938,
               0.42480000853538513, 0.4255000054836273, 0.4207000136375427, 0.4794999957084656, 0.4399999976158142,
               0.44769999384880066, 0.5174999833106995, 0.15629999339580536, 0.44209998846054077, 0.4932999908924103,
               0.4903999865055084]]],
    'chr3': [[94842, 196812807, 1.9900000095367432,
              [0.47519999742507935, 0.42309999465942383, 0.5034999847412109, 0.49779999256134033, 0.4634000062942505,
               0.4875999987125397, 0.4848000109195709, 0.5054000020027161, 0.4823000133037567, 0.4625999927520752,
               0.5281000137329102, 0.550000011920929, 0.46869999170303345]]],
    'chr4': [[29134, 191153264, 1.9600000381469727,
              [0.4104999899864197, 0.5218999981880188, 0.49950000643730164, 0.45730000734329224, 0.4277999997138977,
               0.6115000247955322, 0.46810001134872437, 0.49160000681877136, 0.451200008392334, 0.4607999920845032,
               0.5008000135421753, 0.5677000284194946, 0.4672999978065491, 0.49050000309944153, 0.41690000891685486,
               0.47189998626708984, 0.476500004529953, 0.4625999927520752]]],
    'chr5': [[218403, 180914250, 2.009999990463257,
              [0.6001999735832214, 0.3828999936580658, 0.690200012922287, 0.3474000036716461, 0.3496000111103058,
               0.37929999232292175, 0.634500002861023, 0.6522000193595886, 0.3812999963760376, 0.3706000089645386,
               0.6550000071525574, 0.3221000075340271, 0.3345000088214874, 0.3458000063896179]]]
}


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.min_germline_af = 0.10
        self.max_somatic_af = 0.4
        self.min_nr_SNPs_per_segment = 10
        self.min_segment_length = 10000000
        self.gnomAD_AF_limit = 0.00001
        self.vaf_baseline = 0.48

        self.vcf = ".tests/units/estimate_ctDNA_fraction/sample1_T.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.filter.snv_hard_filter_umi.codon_snvs.sorted.vep_annotated.qci.vcf"  # noqa
        self.germline_vcf = ".tests/units/estimate_ctDNA_fraction/sample1_T.ensembled.vep_annotated.filter.germline.exclude.blacklist.vcf.gz"  # noqa
        self.segments = ".tests/units/estimate_ctDNA_fraction/sample1_T.jumble.pathology_purecn.vcf"
        self.ctDNA_fraction = ".tests/units/estimate_ctDNA_fraction/sample1.ctDNA_fraction.tsv"

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _read_segments(self, test_table, segment_dict):
        for chrom in segment_dict:
            i = 0
            for segment in segment_dict[chrom]:
                try:
                    self.assertEqual(test_table[chrom][i], segment)
                except AssertionError as e:
                    print(f"Failed reading segments. {chrom} {segment} {test_table[chrom][i]}")
                    raise e
                i += 1

    def test_read_segments(self):
        from estimate_ctDNA_fraction import read_segments

        segment_dict = read_segments(self.segments)

        self._read_segments(test_segment_dict, segment_dict)

    def _test_read_germline_vcf(self, test_table, updated_segment_dict):
        for chrom in updated_segment_dict:
            i = 0
            for segment in updated_segment_dict[chrom]:
                try:
                    self.assertEqual(test_table[chrom][i], segment)
                except AssertionError as e:
                    print(f"Failed reading germline vcf. {chrom} {segment} {test_table[chrom][i]}")
                    raise e
                i += 1

    def test_read_germline_vcf(self):
        from estimate_ctDNA_fraction import read_germline_vcf

        updated_segment_dict, germline_dict = read_germline_vcf(self.germline_vcf, test_segment_dict2, self.min_germline_af)
        print(updated_segment_dict)

        self._test_read_germline_vcf(test_updated_dict, updated_segment_dict)

    def test_test_if_signal_in_segment(self):
        from estimate_ctDNA_fraction import test_if_signal_in_segment

        # Test no signal in segment
        signal_bool = test_if_signal_in_segment(test_segment_dict3["chr5"][0][3], 1, 0, self.vaf_baseline)

        test_signal_bool = False

        try:
            self.assertEqual(test_signal_bool, signal_bool)
        except AssertionError as e:
            print(f"Failed testing signal in segment. {test_signal_bool} {signal_bool}")
            raise e

        # Test signal in segment
        signal_bool = test_if_signal_in_segment(test_segment, 1, 0, self.vaf_baseline)

        test_signal_bool = True

        try:
            self.assertEqual(test_signal_bool, signal_bool)
        except AssertionError as e:
            print(f"Failed testing signal in segment. {test_signal_bool} {signal_bool}")
            raise e

    def test_baf_to_tc(self):
        from estimate_ctDNA_fraction import baf_to_tc

        # Test no aneuploidy based on copy number and no BAF signal
        tc, CNA_type = baf_to_tc(0, 2.0, [1.9, 2.0, 2.1], 0)

        test_tc = 0.0
        test_CNA_type = "Del"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

        # Test no aneuploidy based on copy number but with small BAF signal
        tc, CNA_type = baf_to_tc(0.02, 1.95, [1.9, 2.0, 2.1], 0)

        test_tc = 0.07692307692307693
        test_CNA_type = "Del"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

        # Test no aneuploidy in sample based on copy number but with strong BAF signal
        tc, CNA_type = baf_to_tc(0.2, 1.9, [1.9, 2.0, 2.1], 0)

        test_tc = 0.5714285714285715
        test_CNA_type = "Del"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

        # Test aneuploidy in sample based on copy number with deletion in segment
        tc, CNA_type = baf_to_tc(0.2, 0.9, [0.9, 2.0, 3.0], 0)

        test_tc = 0.5714285714285715
        test_CNA_type = "Del"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

        # Test aneuploidy in sample based on copy number with copy neutral LoH in segment
        tc, CNA_type = baf_to_tc(0.2, 1.9, [1.0, 2.0, 3.0], 0)

        test_tc = 0.46136101499423304
        test_CNA_type = "CNLoH"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

        # Test aneuploidy in sample based on copy number with duplication in segment
        tc, CNA_type = baf_to_tc(0.2, 2.9, [1.0, 2.0, 3.0], 0)

        test_tc = 0.7504690431519699
        test_CNA_type = "Dup"

        try:
            self.assertEqual(test_tc, tc)
            self.assertEqual(test_CNA_type, CNA_type)
        except AssertionError as e:
            print(f"Failed calculating tc. {test_tc} {tc}, {test_CNA_type} {CNA_type}")
            raise e

    def test_calculate_cnv_tc(self):
        from estimate_ctDNA_fraction import calculate_cnv_tc

        # No CNA signal found
        tc, seg_list = calculate_cnv_tc(
            test_segment_dict3, self.min_nr_SNPs_per_segment, self.vaf_baseline, self.min_segment_length
        )

        test_tc = 0

        try:
            self.assertEqual(test_tc, tc)
        except AssertionError as e:
            print(f"Failed calculating max tc. {test_tc} {tc}")
            raise e

        # Found a deletion
        tc, seg_list = calculate_cnv_tc(
            test_segment_dict4, self.min_nr_SNPs_per_segment, self.vaf_baseline, self.min_segment_length
        )

        test_tc = 0.36253480233650015

        try:
            self.assertEqual(test_tc, tc)
        except AssertionError as e:
            print(f"Failed calculating max tc. {test_tc} {tc}")
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
            "CHIP_genes": ["", ["DNMT3A", "TET2", "ASXL1"]],
            "Other": ["", []]
        }

        AF, snv_list = read_snv_vcf_and_find_max_af(self.vcf, filter_dict)
        print(AF, snv_list)

        # Updated expected AF based on new filtering logic and substitute VCF
        test_AF = 2.8200000524520874

        try:
            self.assertEqual(test_AF, AF*100)
        except AssertionError as e:
            print(f"Failed reading vcf. {test_AF} {AF*100}")
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
