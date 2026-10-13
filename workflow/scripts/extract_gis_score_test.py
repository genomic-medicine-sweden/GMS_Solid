import pandas as pd
import unittest
import tempfile
import os
import sys
from types import SimpleNamespace

# Add the script directory to sys.path to import it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from extract_gis_score import main, normalise_tc, read_ichorcna, review  # noqa: E402


class TestExtractGisScore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.gis_csv = os.path.join(self.temp_dir.name, "jumble_gis.csv")
        self.output_txt = os.path.join(self.temp_dir.name, "predicted_gis.txt")
        self.purecn = os.path.join(self.temp_dir.name, "purity.txt")
        self.ichorcna = os.path.join(self.temp_dir.name, "params.txt")

        # Create a sample GIS CSV
        df = pd.DataFrame({
            'fraction': [0.1, 0.2, 0.3, 0.4, 0.5],
            'predicted_gis': [10, 20, 30, 40, 50]
        })
        df.to_csv(self.gis_csv, index=False)

    def tearDown(self):
        self.temp_dir.cleanup()

    def run_main(self, tc, pathology_tc="0", purecn="", ichorcna_tf="0", gis_csv=None):
        with open(self.purecn, "w") as f:
            f.write(purecn)
        with open(self.ichorcna, "w") as f:
            f.write(f"Gender:\tfemale\nTumor Fraction:\t{ichorcna_tf}\nPloidy:\t2\n")
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=gis_csv or self.gis_csv, purecn_purity=self.purecn, ichorcna_params=self.ichorcna),
            params=SimpleNamespace(tc=tc, pathology_tc=pathology_tc, hrd_cutoff=42, review_margin=5),
            output=SimpleNamespace(gis_score=self.output_txt),
            wildcards=SimpleNamespace(tc_method="pathology_purecn"),
        )
        main(snakemake)
        return pd.read_csv(self.output_txt, dtype=str, keep_default_na=False).set_index("method")

    def test_header_and_rows(self):
        out = self.run_main(0.3)
        self.assertEqual(list(out.reset_index().columns), ["method", "TC", "predicted_gis", "manual_review"])
        self.assertEqual(list(out.index), ["reported (pathology_purecn)", "pathology", "PureCN", "ichorCNA"])

    def test_extract_exact_match(self):
        out = self.run_main(0.3)
        self.assertEqual(out.loc["reported (pathology_purecn)", "TC"], "0.3")
        self.assertEqual(out.loc["reported (pathology_purecn)", "predicted_gis"], "30")

    def test_extract_nearest_match(self):
        out = self.run_main(0.24)
        self.assertEqual(out.loc["reported (pathology_purecn)", "TC"], "0.24")
        self.assertEqual(out.loc["reported (pathology_purecn)", "predicted_gis"], "20")

    def test_extract_missing_tc(self):
        out = self.run_main(None)
        self.assertEqual(out.loc["reported (pathology_purecn)", "TC"], "-")
        self.assertEqual(out.loc["reported (pathology_purecn)", "predicted_gis"], "-")

    def test_extract_invalid_tc(self):
        out = self.run_main("invalid")
        self.assertEqual(out.loc["reported (pathology_purecn)", "TC"], "invalid")
        self.assertEqual(out.loc["reported (pathology_purecn)", "predicted_gis"], "-")

    def test_extract_empty_df(self):
        empty_csv = os.path.join(self.temp_dir.name, "empty_gis.csv")
        pd.DataFrame(columns=['fraction', 'predicted_gis']).to_csv(empty_csv, index=False)
        out = self.run_main(0.5, pathology_tc="0.4", purecn="0.3", gis_csv=empty_csv)
        self.assertEqual(out.loc["reported (pathology_purecn)", "TC"], "0.5")
        self.assertTrue((out["predicted_gis"] == "-").all())
        self.assertEqual(out.loc["reported (pathology_purecn)", "manual_review"], "no")

    def test_gis_at_each_estimate(self):
        # Pathology in percent, PureCN and ichorCNA as fractions
        out = self.run_main(0.5, pathology_tc="50", purecn="0.3\n", ichorcna_tf="0.2104")
        self.assertEqual(out.loc["pathology", "TC"], "0.5")
        self.assertEqual(out.loc["pathology", "predicted_gis"], "50")
        self.assertEqual(out.loc["PureCN", "predicted_gis"], "30")
        self.assertEqual(out.loc["ichorCNA", "TC"], "0.21")
        self.assertEqual(out.loc["ichorCNA", "predicted_gis"], "20")

    def test_missing_estimates_are_na(self):
        out = self.run_main(0.3, pathology_tc="0", purecn="", ichorcna_tf="0")
        for method in ("pathology", "PureCN", "ichorCNA"):
            self.assertEqual(out.loc[method, "TC"], "-")
            self.assertEqual(out.loc[method, "predicted_gis"], "-")

    def test_manual_review_flag(self):
        # Reported 50 (positive) but PureCN 30 (negative): different HRD calls
        out = self.run_main(0.5, pathology_tc="0.5", purecn="0.3")
        self.assertEqual(out.loc["reported (pathology_purecn)", "manual_review"], "yes")
        # All estimates clearly below the cutoff
        out = self.run_main(0.2, pathology_tc="0.2", purecn="0.1")
        self.assertEqual(out.loc["reported (pathology_purecn)", "manual_review"], "no")
        # The flag is only given on the reported row
        self.assertTrue((out.drop(index="reported (pathology_purecn)")["manual_review"] == "-").all())


class TestHelpers(unittest.TestCase):
    def test_normalise_tc(self):
        self.assertEqual(normalise_tc("0.35"), 0.35)
        self.assertEqual(normalise_tc(35), 0.35)
        for value in (None, "", "NA", "nan", 0, "0", -1):
            self.assertIsNone(normalise_tc(value), value)

    def test_read_ichorcna_both_layouts(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("Sample\tTumor Fraction\tPloidy\nS_T\t0.2604\t2.945\n\n\nS_T\nGender:\tfemale\nTumor Fraction:\t0.2604\n")
        try:
            self.assertEqual(read_ichorcna(f.name), 0.2604)
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
