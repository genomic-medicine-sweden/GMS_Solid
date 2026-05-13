import pandas as pd
import unittest
import tempfile
import os
import sys
from types import SimpleNamespace

# Add the script directory to sys.path to import it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from extract_gis_score import main  # noqa: E402


class TestExtractGisScore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.gis_csv = os.path.join(self.temp_dir.name, "jumble_gis.csv")
        self.output_txt = os.path.join(self.temp_dir.name, "predicted_gis.txt")

        # Create a sample GIS CSV
        df = pd.DataFrame({
            'fraction': [0.1, 0.2, 0.3, 0.4, 0.5],
            'predicted_gis': [10, 20, 30, 40, 50]
        })
        df.to_csv(self.gis_csv, index=False)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_extract_exact_match(self):
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=self.gis_csv),
            params=SimpleNamespace(tc=0.3),
            output=SimpleNamespace(gis_score=self.output_txt)
        )
        main(snakemake)

        with open(self.output_txt, 'r') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "TC,predicted_gis\n")
            self.assertEqual(lines[1], "0.3,30\n")

    def test_extract_nearest_match(self):
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=self.gis_csv),
            params=SimpleNamespace(tc=0.24),
            output=SimpleNamespace(gis_score=self.output_txt)
        )
        main(snakemake)

        with open(self.output_txt, 'r') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "TC,predicted_gis\n")
            self.assertEqual(lines[1], "0.24,20\n")

    def test_extract_missing_tc(self):
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=self.gis_csv),
            params=SimpleNamespace(tc=None),
            output=SimpleNamespace(gis_score=self.output_txt)
        )
        main(snakemake)

        with open(self.output_txt, 'r') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "TC,predicted_gis\n")
            self.assertEqual(lines[1], "NA,NA\n")

    def test_extract_invalid_tc(self):
        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=self.gis_csv),
            params=SimpleNamespace(tc="invalid"),
            output=SimpleNamespace(gis_score=self.output_txt)
        )
        main(snakemake)

        with open(self.output_txt, 'r') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "TC,predicted_gis\n")
            self.assertEqual(lines[1], "invalid,NA\n")

    def test_extract_empty_df(self):
        empty_csv = os.path.join(self.temp_dir.name, "empty_gis.csv")
        pd.DataFrame(columns=['fraction', 'predicted_gis']).to_csv(empty_csv, index=False)

        snakemake = SimpleNamespace(
            input=SimpleNamespace(gis_csv=empty_csv),
            params=SimpleNamespace(tc=0.5),
            output=SimpleNamespace(gis_score=self.output_txt)
        )
        main(snakemake)

        with open(self.output_txt, 'r') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "TC,predicted_gis\n")
            self.assertEqual(lines[1], "0.5,NA\n")


if __name__ == "__main__":
    unittest.main()
