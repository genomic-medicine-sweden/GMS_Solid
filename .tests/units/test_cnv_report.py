import os
import sys
import tempfile
import unittest
from dataclasses import dataclass

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from cnv_report import create_tsv_report  # noqa


class TestGetCaller(unittest.TestCase):
    def test_create_tsv_report(self):
        cnv = tempfile.mkdtemp() + "/tcnv.txt"
        out_additional_only = open(tempfile.mkdtemp() + "/out_additional_only.txt", "w")
        out_tsv_chrom_arms = tempfile.mkdtemp() + "/cnv_chromosome_arms.tsv"
        out_vcf_filename = tempfile.mkdtemp() + "/out_vcf.vcf"
        create_tsv_report(
            [".tests/units/vcf/test.cnv1.vcf"],
            [".tests/units/vcf/test.cnv2.vcf"],
            ".tests/units/vcf/test.deletions.tsv",
            ".tests/units/vcf/test.amplifications.tsv",
            ".tests/integration/reference/chromosome_arm_size.tsv",
            2.0,
            cnv,
            out_additional_only,
            out_tsv_chrom_arms,
            out_vcf_filename,
            1.5,
            0.5,
            0.3,
            1.4,
            2.5,
            1.7,
            2.25,
            0.4,
            0.6,
            0.2,
            0.2,
            0.5,
        )

        @dataclass
        class TestCase:
            name: str
            expected: str

        testcases = [
                TestCase(
                    name="header row",
                    expected=("gene(s)", "chrom", "region", "caller", "freq_in_db", "copy_number")
                ),
                TestCase(
                    name="small deletion",
                    expected=("CDKN2A,CDKN2B", "chr9", "21968207-22008972", "small_deletion", "NA", "-0.28")
                ),
                TestCase(
                    name="variant 1",
                    expected=("FGFR1", "chr8", "34370199-43930232", "cnvkit", "0.01", "8.59")
                ),
                TestCase(
                    name="variant 2",
                    expected=("FGFR1,MYC", "chr8", "35008818-146144253", "gatk", "0.01", "7.01")
                ),
                TestCase(
                    name="variant 3",
                    expected=("MYC", "chr8", "46689525-146144003", "cnvkit", "0.09", "5.06")
                ),
                TestCase(
                    name="small amplification",
                    expected=("MYCN", "chr2", "16968207-17008972", "small_amplification", "NA", "7.29")
                ),
        ]

        testcases2 = [
                TestCase(
                    name="header row",
                    expected=('chrom', 'arm', 'caller', 'type', 'fraction')
                ),
                TestCase(
                    name="baseline 1",
                    expected=("Warning: baseline of GATK CNV might be shifted!", '', '', '', "4.5% on baseline")
                ),
                TestCase(
                    name="baseline 2",
                    expected=("Warning: baseline of CNVkit might be shifted!", '', '', '', "0.0% on baseline")
                ),
                TestCase(
                    name="duplication",
                    expected=('chr8', 'q', 'cnvkit', 'duplication', '98.7%')
                ),
                TestCase(
                    name="loh 1",
                    expected=('chr9', 'p', 'cnvkit', 'loh', '97.5%')
                ),
                TestCase(
                    name="loh 2",
                    expected=('chr9', 'q', 'cnvkit', 'loh', '99.2%')
                ),
        ]

        test_it = iter(testcases)
        with open(cnv) as reader:
            for line in reader:

                actual = tuple(line.rstrip().split("\t"))
                try:
                    case = next(test_it)
                    self.assertEqual(
                        case.expected,
                        actual,
                        "failed test '{}': expected {}, got {}".format(
                            case.name, case.expected, actual
                        ),
                    )
                except ValueError:
                    if case.expected == "ValueError":
                        assert True
                    else:
                        assert False
                except StopIteration:
                    print("More results than test cases")
                    assert False

        test_it = iter(testcases2)
        with open(out_tsv_chrom_arms) as reader:
            for line in reader:

                actual = tuple(line.rstrip().split("\t"))
                try:
                    case = next(test_it)
                    self.assertEqual(
                        case.expected,
                        actual,
                        "failed test '{}': expected {}, got {}".format(
                            case.name, case.expected, actual
                        ),
                    )
                except ValueError:
                    if case.expected == "ValueError":
                        assert True
                    else:
                        assert False
                except StopIteration:
                    print("More results than test cases")
                    assert False
