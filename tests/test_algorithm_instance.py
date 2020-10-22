import os
import sys

import pandas as pd
import pytest

# Add parent directory to search path, so we can import those files
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from algorithm_instance import AlgorithmInstance


def test_runs_without_error():
    example_file = os.path.join(
        os.path.dirname(__file__),
        '../Angola2017_example.xlsx')
    test = AlgorithmInstance(example_file, [2,2,3,3,3,3,3])
    # TODO:
    # test.run_algorithm(nruns=10, burnin=1, record_interval=2, seed=2020)


def test_getting_site_samples():
    expected_sample_df = pd.DataFrame(
        [[
            'BQ17-269_ Day 0', 223.4
        ]],
        columns=['Sample ID', 'Other Data'])

    sample_df = pd.DataFrame(
        [
            ['BQ17-269_ Day 0', 'Benguela', 223.4],
            ['ZL17-019_ Day 0', 'Zaire', 221.3]
        ],
        columns=['Sample ID', 'Site', 'Other Data'])
    sample_df = AlgorithmInstance._get_samples_from_site(sample_df, 'Benguela')

    pd.testing.assert_frame_equal(sample_df, expected_sample_df)


def test_replacing_sample_names():
    expected_sample_df = pd.DataFrame(
        [[
            'Sample_0', 223.4
        ]],
        columns=['Sample ID', 'Other Data'])

    sample_df = pd.DataFrame(
        [[
            'BQ17-269_ Day 0', 223.4
        ]],
        columns=['Sample ID', 'Other Data'])
    AlgorithmInstance._replace_sample_names(sample_df, 'Sample_')

    pd.testing.assert_frame_equal(sample_df, expected_sample_df)

@pytest.mark.xfail(reason='Test not implemented')
def test_calculating_summary_stats():
    # TODO: Implement this!
    assert 'Test not implemented' == ''
