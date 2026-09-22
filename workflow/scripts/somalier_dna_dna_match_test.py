import os
import pandas as pd
import pytest
from somalier_dna_dna_match import main


def test_somalier_dna_dna_match(tmp_path):
    # Create mock Somalier pairs data
    # sample_a and sample_b are suffixed with _T, _N (DNA) or _R (RNA)
    data = {
        '#sample_a': ['S1_T', 'S1_T', 'S1_T', 'S2_T', 'S2_T', 'S3_R'],
        'sample_b': ['S1_T', 'S2_T', 'S3_R', 'S1_T', 'S3_R', 'S1_T'],
        'relatedness': [1.0, 0.8, 0.5, 0.8, 0.4, 0.5],
        'ibs0': [0, 1, 50, 1, 60, 50],
        'ibs2': [100, 80, 40, 80, 30, 40],
        'n': [1000, 1000, 1000, 1000, 1000, 1000]
    }
    df = pd.DataFrame(data)
    input_file = tmp_path / "pairs.tsv"
    df.to_csv(input_file, sep='\t', index=False)

    output_file = tmp_path / "report.tsv"
    match_cutoff = 0.7
    ibs0_cutoff = 2

    # Run the main function
    main(str(input_file), str(output_file), match_cutoff, ibs0_cutoff)

    # Verify results
    assert os.path.exists(output_file)
    res_df = pd.read_csv(output_file, sep='\t')

    # Expected:
    # S1_T best match is S2_T (relatedness 0.8, ibs0 1) -> Match: yes
    # S2_T best match is S1_T (relatedness 0.8, ibs0 1) -> Match: yes
    # S3_R should not be in the results as it is RNA

    assert len(res_df) == 2
    assert 'S1_T' in res_df['Sample'].values
    assert 'S2_T' in res_df['Sample'].values
    assert 'S3_R' not in res_df['Sample'].values

    s1_row = res_df[res_df['Sample'] == 'S1_T'].iloc[0]
    assert s1_row['Best_Match'] == 'S2_T'
    assert s1_row['Relatedness'] == 0.8
    assert s1_row['Match'] == 'yes'

    s2_row = res_df[res_df['Sample'] == 'S2_T'].iloc[0]
    assert s2_row['Best_Match'] == 'S1_T'
    assert s2_row['Relatedness'] == 0.8
    assert s2_row['Match'] == 'yes'


def test_somalier_dna_dna_match_high_ibs0_not_flagged(tmp_path):
    # Relatedness clears match_cutoff, but ibs0 is far from zero, so the AND
    # requirement should keep this from being flagged as a match.
    data = {
        '#sample_a': ['S1_T', 'S2_T'],
        'sample_b': ['S2_T', 'S1_T'],
        'relatedness': [0.8, 0.8],
        'ibs0': [50, 50],
        'ibs2': [80, 80],
        'n': [1000, 1000]
    }
    df = pd.DataFrame(data)
    input_file = tmp_path / "pairs.tsv"
    df.to_csv(input_file, sep='\t', index=False)

    output_file = tmp_path / "report_high_ibs0.tsv"
    main(str(input_file), str(output_file), 0.7, 2)

    assert os.path.exists(output_file)
    res_df = pd.read_csv(output_file, sep='\t')

    assert len(res_df) == 2
    assert (res_df['Match'] == 'no').all()


def test_somalier_dna_dna_match_no_matches(tmp_path):
    # Only self matches or cross-type matches
    data = {
        '#sample_a': ['S1_T', 'S1_T', 'S2_R'],
        'sample_b': ['S1_T', 'S2_R', 'S1_T'],
        'relatedness': [1.0, 0.5, 0.5],
        'ibs0': [0, 50, 50],
        'ibs2': [100, 40, 40],
        'n': [1000, 1000, 1000]
    }
    df = pd.DataFrame(data)
    input_file = tmp_path / "pairs.tsv"
    df.to_csv(input_file, sep='\t', index=False)

    output_file = tmp_path / "report_no_matches.tsv"
    main(str(input_file), str(output_file), 0.7, 2)

    assert os.path.exists(output_file)
    res_df = pd.read_csv(output_file, sep='\t')
    assert len(res_df) == 0


if __name__ == "__main__":
    # This allows running the test script directly if needed
    pytest.main([__file__])
