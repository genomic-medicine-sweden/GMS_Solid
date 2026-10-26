import pandas as pd
import unittest
import tempfile
import os
import sys
from types import SimpleNamespace

# Add the script directory to sys.path to import it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from extract_gis_score import dna_fraction, main, normalise_tc, read_ichorcna, read_purecn, review  # noqa: E402

BEST = "best estimate (PureCN)"


class TestExtractGisScore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.gis_csv = os.path.join(self.temp_dir.name, "jumble_gis.csv")
        self.output_txt = os.path.join(self.temp_dir.name, "predicted_gis.txt")
        self.purecn = os.path.join(self.temp_dir.name, "purecn.csv")
        self.ichorcna = os.path.join(self.temp_dir.name, "params.txt")

        # GIS falls with fraction, one row per 0.05
        fractions = [round(0.05 * i, 2) for i in range(1, 21)]
        pd.DataFrame({"fraction": fractions, "predicted_gis": [round(100 - 100 * f) for f in fractions]}).to_csv(
            self.gis_csv, index=False
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def run_main(self, pathology_tc="0", purecn=(None, None), ichorcna=("0", "2"), gis_csv=None):
        with open(self.purecn, "w") as f:
            f.write('"Sampleid","Purity","Ploidy","Sex"\n')
            if purecn[0] is not None:
                f.write(f'"S_T",{purecn[0]},{purecn[1]},"F"\n')
        with open(self.ichorcna, "w") as f:
            f.write(f"Gender:\tfemale\nTumor Fraction:\t{ichorcna[0]}\nPloidy:\t{ichorcna[1]}\n")
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=gis_csv or self.gis_csv, purecn_csv=self.purecn, ichorcna_params=self.ichorcna),
            params=SimpleNamespace(pathology_tc=pathology_tc, hrd_cutoff=42, review_margin=5, dna_fraction_review_margin=0.1),
            output=SimpleNamespace(gis_score=self.output_txt),
        )
        main(snakemake)
        return pd.read_csv(self.output_txt, dtype=str, keep_default_na=False).set_index("method")

    def test_header_and_rows(self):
        out = self.run_main(purecn=(0.3, 2))
        self.assertEqual(
            list(out.reset_index().columns),
            ["method", "purity", "ploidy", "DNA_fraction", "GIS_ploidy_2", "GIS_DNA_fraction", "manual_review"],
        )
        self.assertEqual(list(out.index), [BEST, "pathology", "PureCN", "ichorCNA"])

    def test_dna_fraction_uses_each_methods_own_ploidy(self):
        out = self.run_main(purecn=(0.4, 3.5), ichorcna=("0.3", "2"))
        # PureCN: 0.4 * 3.5 / (0.4 * 3.5 + 2 * 0.6) = 0.54
        self.assertEqual(out.loc["PureCN", "DNA_fraction"], "0.54")
        self.assertEqual(out.loc["PureCN", "GIS_ploidy_2"], "60")
        self.assertEqual(out.loc["PureCN", "GIS_DNA_fraction"], "45")
        # Ploidy 2: DNA fraction equals purity
        self.assertEqual(out.loc["ichorCNA", "DNA_fraction"], "0.3")
        self.assertEqual(out.loc["ichorCNA", "GIS_ploidy_2"], out.loc["ichorCNA", "GIS_DNA_fraction"])

    def test_best_estimate_is_purecn(self):
        out = self.run_main(pathology_tc="0.6", purecn=(0.4, 3.5), ichorcna=("0.45", "3"))
        for column in ("purity", "ploidy", "DNA_fraction", "GIS_ploidy_2", "GIS_DNA_fraction"):
            self.assertEqual(out.loc[BEST, column], out.loc["PureCN", column])

    def test_pathology_has_no_ploidy(self):
        out = self.run_main(pathology_tc="50", purecn=(0.5, 2))
        self.assertEqual(out.loc["pathology", "purity"], "0.5")
        self.assertEqual(out.loc["pathology", "GIS_ploidy_2"], "50")
        for column in ("ploidy", "DNA_fraction", "GIS_DNA_fraction"):
            self.assertEqual(out.loc["pathology", column], "-")

    def test_missing_estimates(self):
        out = self.run_main(pathology_tc="0", purecn=(None, None), ichorcna=("0", "2"))
        for method in (BEST, "pathology", "PureCN", "ichorCNA"):
            for column in ("purity", "DNA_fraction", "GIS_ploidy_2", "GIS_DNA_fraction"):
                self.assertEqual(out.loc[method, column], "-")
        self.assertEqual(out.loc[BEST, "manual_review"], "no")

    def test_empty_gis_table(self):
        empty_csv = os.path.join(self.temp_dir.name, "empty_gis.csv")
        pd.DataFrame(columns=["fraction", "predicted_gis"]).to_csv(empty_csv, index=False)
        out = self.run_main(pathology_tc="0.4", purecn=(0.3, 2), gis_csv=empty_csv)
        self.assertEqual(out.loc["PureCN", "purity"], "0.3")
        self.assertTrue((out["GIS_ploidy_2"] == "-").all())
        self.assertTrue((out["GIS_DNA_fraction"] == "-").all())

    def test_manual_review_flag(self):
        # Same call, far from the cutoff, DNA fractions agree
        out = self.run_main(pathology_tc="0.2", purecn=(0.2, 2), ichorcna=("0.25", "2"))
        self.assertEqual(out.loc[BEST, "manual_review"], "no")
        # Different HRD calls between pathology (ploidy 2) and PureCN (DNA fraction)
        out = self.run_main(pathology_tc="0.2", purecn=(0.7, 2), ichorcna=("0.7", "2"))
        self.assertEqual(out.loc[BEST, "manual_review"], "yes")
        # Same calls, but PureCN and ichorCNA DNA fractions differ by more than 0.1
        out = self.run_main(pathology_tc="0.2", purecn=(0.2, 2), ichorcna=("0.2", "4"))
        self.assertEqual(out.loc[BEST, "manual_review"], "yes")
        # Only the best estimate row carries the flag
        self.assertTrue((out.drop(index=BEST)["manual_review"] == "-").all())


class TestHelpers(unittest.TestCase):
    def test_normalise_tc(self):
        self.assertEqual(normalise_tc("0.35"), 0.35)
        self.assertEqual(normalise_tc(35), 0.35)
        for value in (None, "", "NA", "nan", 0, "0", -1):
            self.assertIsNone(normalise_tc(value), value)

    def test_dna_fraction(self):
        self.assertEqual(dna_fraction(0.4, 2), 0.4)
        self.assertAlmostEqual(dna_fraction(0.4, 4), 0.5714, places=4)
        self.assertIsNone(dna_fraction(None, 2))
        self.assertIsNone(dna_fraction(0.4, None))

    def test_read_purecn(self):
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            f.write('"Sampleid","Purity","Ploidy","Sex"\n"S_T",0.39,4.43894694488699,"F"\n')
        try:
            purity, ploidy = read_purecn(f.name)
            self.assertEqual(purity, 0.39)
            self.assertAlmostEqual(ploidy, 4.4389, places=4)
        finally:
            os.unlink(f.name)
        self.assertEqual(read_purecn("/nonexistent/purecn.csv"), (None, None))

    def test_read_ichorcna_both_layouts(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("Sample\tTumor Fraction\tPloidy\nS_T\t0.2604\t2.945\n\n\nS_T\nGender:\tfemale\nTumor Fraction:\t0.2604\n"
                    "Ploidy:\t2.945\n")
        try:
            self.assertEqual(read_ichorcna(f.name), (0.2604, 2.945))
        finally:
            os.unlink(f.name)

    def test_review(self):
        # Different HRD calls
        self.assertEqual(review([14, 51, 58], 42, 5), "yes")
        # Same call, but close to the cutoff
        self.assertEqual(review([45, 52], 42, 5), "yes")
        # 42 counts as positive, so 42 and 60 give the same call, but 42 is close to the cutoff
        self.assertEqual(review([42, 60], 42, 5), "yes")
        self.assertEqual(review([36, 30], 42, 5), "no")
        self.assertEqual(review([64, 81, 83], 42, 5), "no")
        self.assertEqual(review([None, None], 42, 5), "no")


if __name__ == "__main__":
    unittest.main()
