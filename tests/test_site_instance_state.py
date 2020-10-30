'''
A series of tests to make sure mcmc.py's functionality is equivalent to mcmc.r's
'''

import os
import sys

import numpy as np
import pandas as pd
import pytest
import scipy.stats as sp_stats

# Add parent directory to search path, so we can import those files
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from algorithm_instance import AlgorithmInstance
from algorithm_site_instance import AlgorithmSiteInstance, SiteInstanceState
from recode_alleles import *
from recrudescence_file_parser import RecrudescenceFileParser
import recrudescence_utils

'''
Full example genotype_RR dataframe in R:

               Sample ID  313_1  313_2  313_3  383_1  383_2  383_3 TA1_1 TA1_2
1        BQ17-269_ Day 0  223.4     NA     NA  103.7  140.3     NA 169.0    NA
2  BQ17-269_ Day Failure  225.5     NA     NA  124.1  162.3     NA 178.0    NA
3        BD17-040_ Day 0  229.3     NA     NA  139.1     NA     NA 175.0    NA
4  BD17-040_ Day Failure  243.6     NA     NA  139.0     NA     NA 175.0    NA
5        BD17-083_ Day 0  221.3     NA     NA  151.6     NA     NA 165.2 180.6
6  BD17-083_ Day Failure  231.5     NA     NA  145.0     NA     NA 162.7 165.6
7        BD17-085_ Day 0  261.7     NA     NA  122.6     NA     NA 162.7    NA
8  BD17-085_ Day Failure  239.6     NA     NA  124.0     NA     NA 162.7    NA
9        BD17-087_ Day 0  261.7     NA     NA  124.0   89.0     NA 162.7    NA
10 BD17-087_ Day Failure  245.3     NA     NA  123.6  140.2     NA 163.0 184.0
11       BD17-090_ Day 0  238.4     NA     NA   87.0  123.9     NA 180.8    NA
12 BD17-090_ Day Failure  219.4  249.7     NA   87.0  123.6     NA 165.8 177.9
   TA1_3 POLYA_1 POLYA_2 POLYA_3 POLYA_4 POLYA_5 PFPK2_1 PFPK2_2 PFPK2_3
1     NA   167.3      NA      NA      NA      NA   149.8   158.0   171.2
2     NA   171.0      NA      NA      NA      NA   162.0   159.0      NA
3     NA   144.0      NA      NA      NA      NA   171.0   168.0      NA
4     NA   164.0      NA      NA      NA      NA   171.0      NA      NA
5     NA   170.2      NA      NA      NA      NA   168.2   177.3   180.3
6     NA   104.8   112.3      NA      NA      NA   168.0   162.0      NA
7     NA   167.1      NA      NA      NA      NA   194.7      NA      NA
8     NA   173.0      NA      NA      NA      NA   189.0      NA      NA
9     NA   167.3      NA      NA      NA      NA   186.0   191.9   194.8
10    NA   155.0      NA      NA      NA      NA   168.2   174.2   177.2
11    NA   151.3      NA      NA      NA      NA   161.9      NA      NA
12    NA   151.3   164.0      NA      NA      NA   161.9   168.1      NA
   PFPK2_4  2490_1  2490_2 TA109_1 TA109_2 TA109_3 TA109_4
1       NA    81.6      NA   164.6   175.0      NA      NA
2       NA    81.8      NA   161.6   172.5      NA      NA
3       NA    82.0      78   149.6   159.9      NA      NA
4       NA    82.0      NA   175.7   187.2      NA      NA
5       NA    81.8      NA   164.2   175.0      NA      NA
6       NA    79.0      NA   149.4   160.0      NA      NA
7       NA    78.6      NA   147.9   158.8      NA      NA
8       NA    88.0      NA   170.0   181.4      NA      NA
9       NA    78.6      NA   148.6   159.9      NA      NA
10      NA    82.0      NA   152.4   163.1      NA      NA
11      NA    81.7      NA   148.8   160.4      NA      NA
12      NA    81.7      NA   163.2   175.3      NA      NA
'''

# NOTE: Makes this reliant on AlgorithmInstance tests passing
example_file = os.path.join(
    os.path.dirname(__file__),
    '../Angola2017_example.xlsx')
genotypedata, additional = RecrudescenceFileParser.parse_file(example_file)
genotypedata_RR = AlgorithmInstance._get_samples_from_site(
    genotypedata, 'Benguela')
additional_neutral = AlgorithmInstance._replace_sample_names(
    AlgorithmInstance._get_samples_from_site(additional, 'Benguela'),
    'Additional_')

expected_maxMOI = 5
locirepeats = np.array([2, 2, 3, 3, 3, 3, 3])
expected_ids = pd.unique(["BQ17-269", "BD17-040", "BD17-083", "BD17-085", "BD17-087", "BD17-090"])
expected_locinames = pd.unique(["313", "383", "TA1", "POLYA", "PFPK2", "2490", "TA109"])
alleles_definitions_RR = AlgorithmSiteInstance._get_allele_definitions(
        genotypedata_RR, additional_neutral, expected_locinames.size, locirepeats)


def test_max_MOI():
    maxMOI = AlgorithmSiteInstance._get_max_MOI(genotypedata_RR)
    assert maxMOI == expected_maxMOI


def test_getting_ids():
    ids = recrudescence_utils.get_sample_ids(genotypedata_RR, 'Day 0')

    np.testing.assert_array_equal(ids, expected_ids)


def test_getting_locinames():
    locinames = pd.unique(genotypedata_RR.columns[1:].str.split("_").str[0])

    np.testing.assert_array_equal(locinames, expected_locinames)


def test_calculate_MOI():
    expected_MOI0 = np.array([3, 2, 3, 2, 3, 2])
    expected_MOIf = np.array([2, 2, 2, 2, 3, 2])

    MOI0, MOIf = SiteInstanceState._calculate_sample_MOI(
        genotypedata_RR, expected_ids, expected_locinames)

    np.testing.assert_array_equal(MOI0, expected_MOI0)
    np.testing.assert_array_equal(MOIf, expected_MOIf)


def test_create_initial_state():
    # TODO: Finish implementing so test will pass
    expected_alleles0_firstCol = np.array([223.4, 229.3, 221.3, 261.7, 261.7, 238.4])
    expected_allelesf_firstCol = np.array([225.5, 243.6, 231.5, 239.6, 245.3, 219.4])
    # NOTE: These values are all 1 less than in the R code since python indexes
    # from 0 (not 1)
    expected_recoded0_firstCol = np.array([4, 11, 1, 3, 3, 7])
    expected_recodedf_firstCol = np.array([5, 2, 6, 8, 12, 0])

    state = SiteInstanceState(expected_ids, expected_locinames, expected_maxMOI, genotypedata_RR, additional_neutral, alleles_definitions_RR)

    np.testing.assert_array_equal(
        state.alleles0[:, 0], expected_alleles0_firstCol)
    np.testing.assert_array_equal(
        state.recoded0[:, 0], expected_recoded0_firstCol)
    np.testing.assert_array_equal(
        state.allelesf[:, 0], expected_allelesf_firstCol)
    np.testing.assert_array_equal(
        state.recodedf[:, 0], expected_recodedf_firstCol)


def test_recode_allele_within_existing_range():
    expected_row_index = 2
    alleles_definitions_RR_subset = np.array([
        [219, 221],
        [221, 223],
        [243, 245],
        [261, 263]
    ])

    row_index = SiteInstanceState._recode_allele(
        alleles_definitions_RR_subset,
        proposed=243.4)

    assert row_index == expected_row_index


def test_recode_allele_proposed_not_in_range():
    alleles_definitions_RR_subset = np.array([
        [219, 221],
        [221, 223],
        [243, 245],
        [261, 263]
    ])

    row_index = SiteInstanceState._recode_allele(
        alleles_definitions_RR_subset,
        proposed=164.3)

    assert np.isnan(row_index)


def test_recode_additional_neutral():
    expected_first_column = np.array([5, 9, 0, 2, 6, 4, 8, 3, 2, 0, 0, 1, 10, 1])

    recoded_additional_neutral = SiteInstanceState.recode_additional_neutral(
        additional_neutral, alleles_definitions_RR, expected_locinames, expected_maxMOI)

    np.testing.assert_array_equal(recoded_additional_neutral.shape, np.array([14, 35]))
    np.testing.assert_array_equal(
        recoded_additional_neutral[:, 0], expected_first_column)


@pytest.mark.xfail(reason='Not sure what to test for')
def test_random_initialization():
    # TODO: How to do expected values for stochastic functions?
    # TODO: Test smaller functions to make sure they're correct first?

    state = SiteInstanceState(expected_ids, expected_locinames, expected_maxMOI, genotypedata_RR, additional_neutral, alleles_definitions_RR)

    state.randomize_initial_assignments(
        expected_ids.size,
        expected_locinames.size,
        expected_maxMOI,
        alleles_definitions_RR,
        2020)


def test_create_dvect():
    alleles_definitions_subset = np.array([
        [217, 263],
        [85, 171],
        [71.5, 203.5], # Largest range (132) should determine dvect size
        [102.5, 180.5],
        [136.5, 202.5],
        [71.5, 89.5],
        [146.5, 197.5],
    ])

    # TODO: What does dvect stand for?
    dvect = SiteInstanceState._get_initial_dvect(alleles_definitions_subset)

    assert dvect.size == 133


def test_get_qq():
    expected_qq = 0.5555556
    hidden0_subset = np.array([
        [0, 1, 1],
        [0, 1, np.nan]
    ])
    hiddenf_subset = np.array([
        [0, 1, np.nan],
        [0, 1, np.nan]
    ])

    qq = SiteInstanceState._get_initial_qq(hidden0_subset, hiddenf_subset)
    np.testing.assert_almost_equal(qq, expected_qq, decimal=5)


def test_correction_factor():
    expected_fist_row_col = np.array([0, 2, 24, 42, 4, 6, 12, 18, 20, 30, 2, 10, 26])

    # TODO: Confirm what the right alleles_definitions_RR shape is?
    alleles_definitions_RR_stub = [np.zeros((13, 2)), np.zeros((13, 2))]
    alleles_definitions_RR_stub[0][:, :] = np.array([
        [219, 221],
        [221, 223],
        [243, 245],
        [261, 263],
        [223, 225],
        [225, 227],
        [231, 233],
        [237, 239],
        [239, 241],
        [249, 251],
        [217, 219],
        [229, 231],
        [245, 247],
    ])

    correction_distance_matrix = SiteInstanceState._get_correction_distances(alleles_definitions_RR_stub)

    assert len(correction_distance_matrix) == 2
    np.testing.assert_array_equal(
        correction_distance_matrix[0].shape, np.array([13, 13]))
    np.testing.assert_array_equal(
        correction_distance_matrix[0][:, 0], expected_fist_row_col)
    np.testing.assert_array_equal(
        correction_distance_matrix[0][0, :], expected_fist_row_col)


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_find_allele_modes():
    maxMOI = 5
    nids = 6
    nloci = 7

    # TODO: Stubbed data
    alleles0 = 100 * np.random.random_sample((nids, maxMOI * nloci))
    allelesf = 100 * np.random.random_sample((nids, maxMOI * nloci))

    state_alleles0 = np.full_like(np.empty((nids, maxMOI * nloci, 12)), np.nan)
    state_allelesf = np.full_like(np.empty((nids, maxMOI * nloci, 12)), np.nan)

    for i in range(12):
        state_alleles0[:, :, i] = alleles0
        state_allelesf[:, :, i] = allelesf

    ## find mode of hidden alleles
    # TODO: Why was this an array of strings?
    modealleles = np.zeros((2 * nids, maxMOI * nloci))
    for i in range(nids):
        for j in range(nloci):
            modealleles[2 * i, j * maxMOI : (j + 1) * maxMOI] = sp_stats.mode(
                state_alleles0[i, j * maxMOI : (j + 1) * maxMOI, :], axis=1
            )[0].ravel()

            modealleles[2 * i + 1, j * maxMOI : (j + 1) * maxMOI] = sp_stats.mode(
                state_allelesf[i, j * maxMOI : (j + 1) * maxMOI, :], axis=1
            )[0].ravel()


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_summary_stats_output():
    maxMOI = 5
    nids = 6
    nloci = 7

    ## TODO: Stubbed data
    locinames = pd.unique(["X313", "X383", "TA1", "POLYA", "PFPK2", "X2490", "TA109"])
    state_parameters = np.random.random_sample((2 + 2 * nloci, 25))
    jobname = "TEST"

    # summary statistics of parameters
    # pd.DataFrame(state_parameters).to_csv(f"{jobname}_state_parameters.csv")

    summary_statisticsmatrix = np.concatenate(
        (
            np.mean(state_parameters, axis=1).reshape(-1, 1),
            np.quantile(state_parameters, (0.25, 0.75), axis=1).T,
        ),
        axis=1,
    )
    summary_statisticsmatrix = np.concatenate(
        (
            summary_statisticsmatrix,
            np.append(
                np.quantile(state_parameters[2 + nloci :, :], (0.25, 0.75)),
                np.mean(state_parameters[2 + nloci :, :]),
            ).reshape(1, -1),
        )
    )
    summary_statisticsmatrix = np.array(
        [
            f"{summary_statisticsmatrix[i,0]:.2f} ({summary_statisticsmatrix[i,1]:.2f}, {summary_statisticsmatrix[i,2]:.2f})"
            for i in range(summary_statisticsmatrix.shape[0])
        ]
    )
    summary_statisticsmatrix_df = pd.DataFrame(
        summary_statisticsmatrix,
        index=["q", "d", *locinames.tolist(), *locinames.tolist(), "Mean diversity"],
    )
    # summary_statisticsmatrix_df.to_csv(f"{jobname}_summarystatistics.csv")
