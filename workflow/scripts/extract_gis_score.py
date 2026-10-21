import csv
import math
import os
import sys

import pandas as pd

HEADER = ["method", "TC", "predicted_gis", "manual_review"]


def normalise_tc(value):
    """Tumor content as a fraction, or None when missing, invalid or 0. Values above 1 are percentages."""
    try:
        tc = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(tc) or tc <= 0:
        return None
    return round(tc / 100, 2) if tc > 1 else tc


def read_purecn(path):
    with open(path) as f:
        return normalise_tc(f.read().strip())


def read_ichorcna(path):
    """Tumor fraction from an ichorCNA params.txt (the 'Tumor Fraction:' line)."""
    with open(path) as f:
        for line in f:
            if line.startswith("Tumor Fraction:"):
                return normalise_tc(line.split(":", 1)[1].strip())
    return None


def gis_at(df, tc):
    """Predicted GIS at the fraction closest to tc, or None without a TC or GIS table."""
    if tc is None or df.empty:
        return None
    return df.loc[(df["fraction"] - tc).abs().idxmin(), "predicted_gis"]


def review(gis_values, cutoff, margin):
    """'yes' when the GIS values would give different HRD calls or any is close to the cutoff, else 'no'."""
    gis_values = [g for g in gis_values if g is not None]
    if not gis_values:
        return "no"
    different_calls = min(gis_values) < cutoff <= max(gis_values)
    close_to_cutoff = any(abs(g - cutoff) <= margin for g in gis_values)
    return "yes" if different_calls or close_to_cutoff else "no"


def fmt(value):
    # "-" rather than NA: the general report reads NA as a missing number and shows NaN
    return "-" if value is None else value


def fmt_tc(tc):
    return "-" if tc is None else round(tc, 2)


def main(snakemake):
    gis_csv = snakemake.input.gis_csv
    tc = snakemake.params.tc
    output = snakemake.output.gis_score

    if not os.path.exists(gis_csv):
        print(f"Error: GIS CSV file {gis_csv} does not exist.", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(gis_csv)
    except Exception as e:
        print(f"Error: Could not read GIS CSV file {gis_csv}: {e}", file=sys.stderr)
        sys.exit(1)

    if not df.empty and ("fraction" not in df.columns or "predicted_gis" not in df.columns):
        print(f"Error: GIS CSV file {gis_csv} is missing required columns 'fraction' or 'predicted_gis'.", file=sys.stderr)
        sys.exit(1)
    if df.empty:
        print(f"Warning: GIS CSV file {gis_csv} is empty. Writing '-' GIS values.", file=sys.stderr)

    # The TC used for reporting (from get_tc). It can be empty if not identified, or not a number.
    if tc is None or tc == "" or (isinstance(tc, float) and pd.isna(tc)):
        print("Warning: Tumor content (TC) is not provided or could not be identified. Writing '-'.", file=sys.stderr)
        reported_tc, reported_label = None, "-"
    else:
        try:
            reported_tc = float(tc)
            reported_label = fmt_tc(reported_tc)
        except ValueError:
            print(f"Error: Tumor content '{tc}' is not a valid number. Writing '-' GIS.", file=sys.stderr)
            reported_tc, reported_label = None, tc
    reported_gis = gis_at(df, reported_tc)

    # GIS at every available TC estimate, so differences between them are visible
    estimates = [
        ("pathology", normalise_tc(snakemake.params.pathology_tc)),
        ("PureCN", read_purecn(snakemake.input.purecn_purity)),
        ("ichorCNA", read_ichorcna(snakemake.input.ichorcna_params)),
    ]
    estimate_gis = [(name, est_tc, gis_at(df, est_tc)) for name, est_tc in estimates]

    flag = review(
        [reported_gis] + [g for _, _, g in estimate_gis],
        int(snakemake.params.hrd_cutoff),
        int(snakemake.params.review_margin),
    )

    with open(output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(HEADER)
        writer.writerow([f"reported ({snakemake.wildcards.tc_method})", reported_label, fmt(reported_gis), flag])
        for name, est_tc, gis in estimate_gis:
            writer.writerow([name, fmt_tc(est_tc), fmt(gis), "-"])


if __name__ == "__main__":
    main(snakemake)
