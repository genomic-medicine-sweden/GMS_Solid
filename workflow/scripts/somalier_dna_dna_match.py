import sys
import pandas as pd
import os
import argparse


def get_sample_type(sample_name):
    """
    Identifies sample type based on suffix:
    _T or _N -> DNA
    _R       -> RNA
    """
    if sample_name.endswith('_T') or sample_name.endswith('_N'):
        return 'DNA'
    elif sample_name.endswith('_R'):
        return 'RNA'
    return 'Unknown'


def main(input_file, output_file, match_cutoff, ibs0_cutoff):
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)

    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        sys.exit(1)

    if df.empty:
        print(f"Warning: Input file {input_file} is empty.")
        # Create empty output with headers
        pd.DataFrame(
            columns=['Sample', 'Best_Match', 'Relatedness', 'ibs0', 'Variants_compared', 'Match']
        ).to_csv(output_file, sep='\t', index=False)
        return

    # Column names might have leading #
    rename_dict = {col: col.lstrip('#') for col in df.columns}
    df = df.rename(columns=rename_dict)

    # Required columns
    required = ['sample_a', 'sample_b', 'relatedness', 'ibs0', 'ibs2', 'n']
    for col in required:
        if col not in df.columns:
            print(f"Error: Missing required column {col}")
            sys.exit(1)

    # Get unique samples
    all_samples = pd.concat([df['sample_a'], df['sample_b']]).unique()

    results = []
    for sample in all_samples:
        sample_type = get_sample_type(sample)

        # Only process DNA samples
        if sample_type != 'DNA':
            continue

        # Matches involving this sample
        matches = df[(df['sample_a'] == sample) | (df['sample_b'] == sample)].copy()

        if matches.empty:
            continue

        # Identify the other sample in the pair
        matches['other_sample'] = matches.apply(
            lambda row: row['sample_b'] if row['sample_a'] == sample else row['sample_a'],
            axis=1
        )

        # Filter for other DNA samples
        filtered_matches = matches[matches['other_sample'].apply(lambda s: get_sample_type(s) == 'DNA')]

        # Exclude self-matches (relatedness should be 1.0, but we want the best *other* match)
        # However, Somalier usually pairs different samples. If same sample is present, exclude it.
        filtered_matches = filtered_matches[filtered_matches['other_sample'] != sample]

        if filtered_matches.empty:
            continue

        # Sort by relatedness descending to find the best match
        best_match_row = filtered_matches.sort_values(by='relatedness', ascending=False).iloc[0]

        relatedness = best_match_row['relatedness']
        ibs0 = best_match_row['ibs0']
        # A pair is only flagged as a match if BOTH conditions hold: relatedness clears
        # match_cutoff AND ibs0 is near zero - ibs0 near 0 (few/no sites with opposite
        # homozygous genotypes) is the more specific signature of a real identity issue
        # (swap or cross-contamination), and requiring both avoids flagging samples that
        # simply happen to be each other's least-bad match without real signal behind it.
        match_status = "yes" if (relatedness >= match_cutoff and ibs0 <= ibs0_cutoff) else "no"

        results.append({
            'Sample': sample,
            'Best_Match': best_match_row['other_sample'],
            'Relatedness': relatedness,
            'ibs0': ibs0,
            'Variants_compared': best_match_row['n'],
            'Match': match_status,
        })

    res_df = pd.DataFrame(results, columns=['Sample', 'Best_Match', 'Relatedness', 'ibs0', 'Variants_compared', 'Match'])
    if not res_df.empty:
        # Sort results by sample name for consistency
        res_df = res_df.sort_values(by='Sample')
    res_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        main(
            snakemake.input.pairs,
            snakemake.output.report,
            float(snakemake.params.match_cutoff),
            float(snakemake.params.ibs0_cutoff),
        )
    else:
        parser = argparse.ArgumentParser(description="Find best DNA-DNA matches from Somalier pairs.")
        parser.add_argument("input", help="Input Somalier pairs TSV file")
        parser.add_argument("output", help="Output best match TSV file")
        parser.add_argument("--match-cutoff", type=float, default=0.6, help="Relatedness cutoff for match (default: 0.6)")
        parser.add_argument(
            "--ibs0-cutoff", type=float, default=2,
            help="Max ibs0 required (together with --match-cutoff) to flag a match (default: 2)",
        )

        args = parser.parse_args()
        main(args.input, args.output, args.match_cutoff, args.ibs0_cutoff)
