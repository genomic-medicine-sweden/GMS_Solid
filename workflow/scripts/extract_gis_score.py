import pandas as pd
import sys
import os


def main(snakemake):
    gis_csv = snakemake.input.gis_csv
    tc = snakemake.params.tc
    output = snakemake.output.gis_score

    # tc can be an empty string if not identified in get_tc
    if tc is None or tc == "" or (isinstance(tc, float) and pd.isna(tc)):
        print("Warning: Tumor content (TC) is not provided or could not be identified. Writing 'NA' to output.", file=sys.stderr)
        with open(output, "w") as f:
            f.write("TC,predicted_gis\n")
            f.write("NA,NA\n")
        return

    try:
        tc_val = float(tc)
    except ValueError:
        print(f"Error: Tumor content '{tc}' is not a valid number. Writing 'NA' to output.", file=sys.stderr)
        with open(output, "w") as f:
            f.write("TC,predicted_gis\n")
            f.write(f"{tc},NA\n")
        return

    if not os.path.exists(gis_csv):
        print(f"Error: GIS CSV file {gis_csv} does not exist.", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(gis_csv)
    except Exception as e:
        print(f"Error: Could not read GIS CSV file {gis_csv}: {e}", file=sys.stderr)
        sys.exit(1)

    if 'fraction' not in df.columns or 'predicted_gis' not in df.columns:
        print(f"Error: GIS CSV file {gis_csv} is missing required columns 'fraction' or 'predicted_gis'.", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print(f"Warning: GIS CSV file {gis_csv} is empty. Writing 'NA' to output.", file=sys.stderr)
        with open(output, "w") as f:
            f.write("TC,predicted_gis\n")
            f.write(f"{tc_val},NA\n")
        return

    # Find the row with the fraction closest to tc_val
    idx = (df['fraction'] - tc_val).abs().idxmin()
    predicted_gis = df.loc[idx, 'predicted_gis']

    with open(output, "w") as f:
        f.write("TC,predicted_gis\n")
        f.write(f"{tc_val},{predicted_gis}\n")


if __name__ == "__main__":
    main(snakemake)
