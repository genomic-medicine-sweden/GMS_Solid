import os
import sys
import tempfile
import unittest

import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ichorcna_purity import main  # noqa: E402

PARAMS_TXT = (
    "Sample\tTumor Fraction\tPloidy\n"
    "S1_T\t0.1234\t2.01\n"
    "\n"
    "\n"
    "S1_T\n"
    "Gender:\tfemale\n"
    "Tumor Fraction:\t0.1234\n"
    "Ploidy:\t2.01\n"
    "Subclone Fraction:\tNA\n"
    "Fraction Genome Subclonal:\t0\n"
    "Fraction CNA Subclonal:\t0\n"
    "Coverage:\tNA\n"
    "ChrY coverage fraction:\t0.0012\n"
    "ChrX median log ratio:\t0.0213\n"
    "Student's t mean: 0.01, -0.5, 0.3\n"
    "Student's t precision: 1500, 1500, 1500\n"
    "Gamma Rate Init:\t1e+05\n"
    "GC-Map correction MAD:\t0.0871\n"
    "\n"
    "\n"
    "init\tn_est\tphi_est\tBIC\tFrac_genome_subclonal\tFrac_CNA_subclonal\tloglik\n"
    "n0.5-p2\t0.88\t2.01\tNA\t0\t0\t-120\n"
    "n0.5-p3\t0.6\t3.1\tNA\t0.9\t0.8\t-100\n"
    "n0.5-p4\t0.75\t3.9\tNA\t0.1\t0.2\t-130\n"
)


class TestIchorcnaPurity(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.params_txt = os.path.join(self.temp_dir.name, "S1_T.params.txt")
        self.output = os.path.join(self.temp_dir.name, "S1_T.ichorcna_purity.tsv")

    def tearDown(self):
        self.temp_dir.cleanup()

    def run_main(self, content):
        with open(self.params_txt, "w") as f:
            f.write(content)
        main(self.params_txt, self.output)
        return pd.read_csv(self.output, sep="\t", keep_default_na=False)

    def test_selected_solution(self):
        df = self.run_main(PARAMS_TXT)
        self.assertEqual(
            list(df.columns),
            ["Sample", "Tumor Fraction", "Ploidy", "Subclone Fraction", "Predicted gender", "GC-Map correction MAD"],
        )
        self.assertEqual(len(df), 1)
        row = df.iloc[0]
        self.assertEqual(row["Sample"], "S1_T")
        self.assertEqual(row["Tumor Fraction"], 0.1234)
        self.assertEqual(row["Ploidy"], 2.01)
        self.assertEqual(row["Subclone Fraction"], "NA")
        self.assertEqual(row["Predicted gender"], "female")
        self.assertEqual(row["GC-Map correction MAD"], 0.0871)

    def test_subclone_fraction_reported(self):
        df = self.run_main(PARAMS_TXT.replace("Subclone Fraction:\tNA", "Subclone Fraction:\t0.3512"))
        self.assertEqual(df.loc[0, "Subclone Fraction"], 0.3512)


if __name__ == "__main__":
    unittest.main()
