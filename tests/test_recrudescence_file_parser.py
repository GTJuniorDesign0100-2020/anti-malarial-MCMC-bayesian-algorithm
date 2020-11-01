import os
import sys

import numpy as np
import pandas as pd

# Add parent directory to search path, so we can import those files
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from api.recrudescence_file_parser import RecrudescenceFileParser


def test_angola_example_parser_output_correct():
    expected_genotype_shape = (76, 26)
    expected_additional_shape = (25, 26)

    expected_columns = np.array([
        'Sample ID', 'Site',
        '313_1', '313_2', '313_3',
        '383_1', '383_2', '383_3',
        'TA1_1', 'TA1_2', 'TA1_3',
        'POLYA_1', 'POLYA_2', 'POLYA_3', 'POLYA_4', 'POLYA_5',
        'PFPK2_1', 'PFPK2_2', 'PFPK2_3', 'PFPK2_4',
        '2490_1', '2490_2',
        'TA109_1', 'TA109_2', 'TA109_3', 'TA109_4'
    ])
    # TODO: Find better way than just testing the 1st row?
    expected_genotype_row_0 = pd.DataFrame(
        [[
            'BQ17-269 Day 0', 'Benguela',
            223.4, np.nan, np.nan,
            103.7, 140.3, np.nan,
            169, np.nan, np.nan,
            167.3, np.nan, np.nan, np.nan, np.nan,
            149.8, 158, 171.2, np.nan,
            81.6, np.nan,
            164.6, 175, np.nan, np.nan
        ]],
        columns=expected_columns)
    expected_additional_row_0 = pd.DataFrame(
        [[
            'BD17-047 Day 0', 'Benguela',
            226, np.nan, np.nan,
            124, 171, np.nan,
            160, 172, np.nan,
            154, 151, np.nan, np.nan, np.nan,
            162, 165, 174, 177,
            82, 78,
            164, 160, np.nan, np.nan
        ]],
        columns=expected_columns)

    example_file = os.path.join(
        os.path.dirname(__file__),
        '../Angola2017_example.xlsx')
    genotype, additional = RecrudescenceFileParser.parse_file(example_file)

    assert genotype.shape == expected_genotype_shape
    assert additional.shape == expected_additional_shape

    np.testing.assert_array_equal(genotype.columns.to_numpy(), expected_columns)
    np.testing.assert_array_equal(additional.columns.to_numpy(), expected_columns)

    pd.testing.assert_series_equal(genotype.iloc[0], expected_genotype_row_0.iloc[0])
    pd.testing.assert_series_equal(additional.iloc[0], expected_additional_row_0.iloc[0])
