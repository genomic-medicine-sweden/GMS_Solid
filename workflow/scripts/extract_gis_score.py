import csv
import math
import os
import sys

import pandas as pd

HEADER = ["method", "purity", "ploidy", "DNA_fraction", "GIS_ploidy_2", "GIS_DNA_fraction", "manual_review"]


def normalise_tc(value):
    """Tumor content as a fraction, or None when missing, invalid or 0. Values above 1 are percentages."""
    try:
        tc = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(tc) or tc <= 0:
        return None
    return round(tc / 100, 2) if tc > 1 else tc


def to_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def read_purecn(path):
    """Purity and ploidy from the same PureCN fit (its csv), or None for each when missing."""
    try:
        row = pd.read_csv(path).iloc[0]
    except (FileNotFoundError, pd.errors.EmptyDataError, IndexError):
        return None, None
    ploidy = to_float(row.get("Ploidy"))
    return normalise_tc(row.get("Purity")), ploidy if ploidy and ploidy > 0 else None


def read_ichorcna(path):
    """Tumor fraction and ploidy from an ichorCNA params.txt ('Tumor Fraction:' and 'Ploidy:' lines)."""
    values = {}
    with open(path) as f:
        for line in f:
            for key in ("Tumor Fraction", "Ploidy"):
                if line.startswith(f"{key}:") and key not in values:
                    values[key] = line.split(":", 1)[1].strip()
    ploidy = to_float(values.get("Ploidy"))
    return normalise_tc(values.get("Tumor Fraction")), ploidy if ploidy and ploidy > 0 else None


def dna_fraction(purity, ploidy):
    """Tumor DNA fraction, the fraction Jumble's GIS table is indexed by, from tumor cell purity and ploidy."""
    if purity is None or ploidy is None:
        return None
    return purity * ploidy / (purity * ploidy + 2 * (1 - purity))


def gis_at(df, fraction):
    """Predicted GIS at the table fraction closest to fraction, or None without a fraction or GIS table."""
    if fraction is None or df.empty:
        return None
    return df.loc[(df["fraction"] - fraction).abs().idxmin(), "predicted_gis"]


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


def fmt_number(value):
    return "-" if value is None else round(value, 2)


def main(snakemake):
    gis_csv = snakemake.input.gis_csv
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

    # Each method's purity and ploidy come from the same fit and are only combined with each other
    methods = [
        ("pathology", normalise_tc(snakemake.params.pathology_tc), None),
        ("PureCN", *read_purecn(snakemake.input.purecn_csv)),
        ("ichorCNA", *read_ichorcna(snakemake.input.ichorcna_params)),
    ]
    rows = []
    for name, purity, ploidy in methods:
        fraction = dna_fraction(purity, ploidy)
        rows.append([name, purity, ploidy, fraction, gis_at(df, purity), gis_at(df, fraction)])

    # Best estimate: PureCN's own purity and ploidy
    best = ["best estimate (PureCN)"] + rows[1][1:]

    pathology, purecn, ichorcna = rows
    dna_fractions_differ = (
        purecn[3] is not None
        and ichorcna[3] is not None
        and abs(purecn[3] - ichorcna[3]) > float(snakemake.params.dna_fraction_review_margin)
    )
    gis_differ = review(
        [best[5], ichorcna[5], pathology[4]], int(snakemake.params.hrd_cutoff), int(snakemake.params.review_margin)
    )
    flag = "yes" if gis_differ == "yes" or dna_fractions_differ else "no"

    with open(output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(HEADER)
        for row, row_flag in [(best, flag)] + [(row, "-") for row in rows]:
            name, purity, ploidy, fraction, gis_ploidy_2, gis_dna = row
            writer.writerow(
                [name, fmt_number(purity), fmt_number(ploidy), fmt_number(fraction), fmt(gis_ploidy_2), fmt(gis_dna), row_flag]
            )


if __name__ == "__main__":
    main(snakemake)
