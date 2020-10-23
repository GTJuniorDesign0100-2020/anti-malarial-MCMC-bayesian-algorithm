"""
A series of tests to make sure mcmc.py's functionality is equivalent to mcmc.r's

TODO: Use an actual testing framework?
TODO: Double-check axes of EVERYTHING (I'm worried I might've mixed up running function on column/row axes with numpy functions)
"""

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

"""
Full example genotype_RR dataframe in R:

               Sample.ID X313_1 X313_2 X313_3 X383_1 X383_2 X383_3 TA1_1 TA1_2
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
   PFPK2_4 X2490_1 X2490_2 TA109_1 TA109_2 TA109_3 TA109_4
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
"""

# Hardcode this in for now
genotypedata_RR = pd.DataFrame(
    {
        "Sample ID": [
            "BQ17-269_ Day 0",
            "BQ17-269_ Day Failure",
            "BD17-040_ Day 0",
            "BD17-040_ Day Failure",
            "BD17-083_ Day 0",
            "BD17-083_ Day Failure",
            "BD17-085_ Day 0",
            "BD17-085_ Day Failure",
            "BD17-087_ Day 0",
            "BD17-087_ Day Failure",
            "BD17-090_ Day 0",
            "BD17-090_ Day Failure",
        ],
        "313_1": [
            223.4,
            225.5,
            229.3,
            243.6,
            221.3,
            231.5,
            261.7,
            239.6,
            261.7,
            245.3,
            238.4,
            219.4,
        ],
        "313_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            249.7,
        ],
        "313_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "383_1": [
            103.7,
            124.1,
            139.1,
            139.0,
            151.6,
            145.0,
            122.6,
            124.0,
            124.0,
            123.6,
            87.0,
            87.0,
        ],
        "383_2": [
            140.3,
            162.3,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            89.0,
            140.2,
            123.9,
            123.6,
        ],
        "383_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA1_1": [
            169.0,
            178.0,
            175.0,
            175.0,
            165.2,
            162.7,
            162.7,
            162.7,
            162.7,
            163.0,
            180.8,
            165.8,
        ],
        "TA1_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            180.6,
            165.6,
            np.nan,
            np.nan,
            np.nan,
            184.0,
            np.nan,
            177.9,
        ],
        "TA1_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_1": [
            167.3,
            171.0,
            144.0,
            164.0,
            170.2,
            104.8,
            167.1,
            173.0,
            167.3,
            155.0,
            151.3,
            151.3,
        ],
        "POLYA_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            112.3,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            164.0,
        ],
        "POLYA_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_5": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "PFPK2_1": [
            149.8,
            162.0,
            171.0,
            171.0,
            168.2,
            168.0,
            194.7,
            189.0,
            186.0,
            168.2,
            161.9,
            161.9,
        ],
        "PFPK2_2": [
            158.0,
            159.0,
            168.0,
            np.nan,
            177.3,
            162.0,
            np.nan,
            np.nan,
            191.9,
            174.2,
            np.nan,
            168.1,
        ],
        "PFPK2_3": [
            171.2,
            np.nan,
            np.nan,
            np.nan,
            180.3,
            np.nan,
            np.nan,
            np.nan,
            194.8,
            177.2,
            np.nan,
            np.nan,
        ],
        "PFPK2_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "2490_1": [
            81.6,
            81.8,
            82.0,
            82.0,
            81.8,
            79.0,
            78.6,
            88.0,
            78.6,
            82.0,
            81.7,
            81.7,
        ],
        "2490_2": [
            np.nan,
            np.nan,
            78,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA109_1": [
            164.6,
            161.6,
            149.6,
            175.7,
            164.2,
            149.4,
            147.9,
            170.0,
            148.6,
            152.4,
            148.8,
            163.2,
        ],
        "TA109_2": [
            175.0,
            172.5,
            159.9,
            187.2,
            175.0,
            160.0,
            158.8,
            181.4,
            159.9,
            163.1,
            160.4,
            175.3,
        ],
        "TA109_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA109_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
    }
)

# TODO: Avoid this, since it makes the test dependent on untested functions
example_file = os.path.join(
    os.path.dirname(__file__),
    '../Angola2017_example.xlsx')
ignore, additional = RecrudescenceFileParser.parse_file(example_file)
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
    expected_recoded0_firstCol = np.array([5, 12, 2, 4, 4, 8])
    expected_recodedf_firstCol = np.array([6, 3, 7, 9, 13, 1])

    state = SiteInstanceState(expected_ids, expected_locinames, expected_maxMOI, genotypedata_RR, additional_neutral, alleles_definitions_RR)

    np.testing.assert_array_equal(
        state.alleles0[:, 0], expected_alleles0_firstCol)
    np.testing.assert_array_equal(
        state.recoded0[:, 0], expected_recoded0_firstCol)
    np.testing.assert_array_equal(
        state.allelesf[:, 0], expected_allelesf_firstCol)
    np.testing.assert_array_equal(
        state.recodedf[:, 0], expected_recodedf_firstCol)


def test_recode_additional_neutral():
    expected_first_column = np.array([6, 10, 1, 3, 7, 5, 9, 4, 3, 1, 0, 2, 11, 2])

    recoded_additional_neutral = SiteInstanceState.recode_additional_neutral(
        additional_neutral, alleles_definitions_RR, expected_locinames, expected_maxMOI)

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


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_create_dvect():
    # TODO: This is a stub, not the actual data
    nloci = 7
    alleles_definitions_RR = np.array(
        [
            [217, 263],
            [85, 171],
            [71.5, 203.5],
            [102.5, 180.5],
            [136.5, 202.5],
            [71.5, 89.5],
            [146.5, 197.5],
        ]
    )

    # TODO: What does dvect stand for?
    dvect = np.zeros(
        1
        + int(
            np.rint(
                # Get the range (max-min) of the first "nloci" rows, then the max of all those
                np.ptp(alleles_definitions_RR[0:nloci], axis=1).max()
            )
        )
    )
    dvect[1] = 0.75
    dvect[2] = 0.2
    dvect[3] = 0.05

    assert dvect.size == 133, f"Dvect size {dvect.size} (expected {133})"


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_correction_factor():
    expected_fist_row_col = np.array([0, 2, 24, 42, 4, 6, 12, 18, 20, 30, 2, 10, 26])

    nloci = 7
    # TODO: Confirm what the right alleles_definitions_RR shape is?
    alleles_definitions_RR = np.zeros((nloci, 13, 2))
    alleles_definitions_RR[0, :, :] = np.array(
        [
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
        ]
    )

    correction_distance_matrix = np.zeros(
        (nloci, alleles_definitions_RR.shape[1], alleles_definitions_RR.shape[1])
    )  # for each locus, matrix of distances between each allele
    for i in range(nloci):
        # Wrap mean call in "array" so we get a 2D array we can transpose (getting us a grid of distances, not just a 1D vector)
        distances = np.array([np.mean(alleles_definitions_RR[i], axis=1)])
        distance_combinations = np.abs(distances.T - distances)
        correction_distance_matrix[i] = distance_combinations

    assert np.array_equal(
        correction_distance_matrix.shape, np.array([7, 13, 13])
    ), f"{correction_distance_matrix.shape} (expected {np.array([7, 13, 13])})"
    assert np.array_equal(
        correction_distance_matrix[0, :, 0], expected_fist_row_col
    ), f"Row {correction_distance_matrix[0, :, 0]} (expected {expected_fist_row_col})"
    assert np.array_equal(
        correction_distance_matrix[0, 0, :], expected_fist_row_col
    ), f"Column {correction_distance_matrix[0, 0, :]} (expected {expected_fist_row_col})"


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_run_mcmc_new_proposal():
    # TODO: Find what the actual expected ratio is for the stubbed inputs?
    expected_likelihood_ratio = np.zeros(6)

    maxMOI = 5
    nids = 6
    nloci = 7

    # TODO: Stubbed data
    frequencies_RR = np.random.random_sample((3, 7, 19))
    frequencies_RR[0, :, 0] = np.array([13, 16, 11, 13, 19, 5, 13])

    alldistance = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), 1)
    allrecrf = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), 3).astype(int)
    hidden0 = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    hiddenf = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    classification = np.repeat(0, nids)

    dvect = np.ones(133)  # TODO: What to do when dvect has a 0 value?
    # TODO: What size is the correction matrix??? Appears to be different for each locus??? (set to max value in frequencies_rr)
    correction_distance_matrix = np.zeros((nloci, 19, 19)).astype(int)

    # propose new classification
    likelihoodratio = np.zeros(nids)
    # TODO: Finish vectorizing this
    for x in range(nids):
        # id mean for what?
        id_means = np.zeros(nloci)
        for y in range(nloci):
            id_means[y] = np.nanmean(
                dvect[np.round(alldistance[x, y, :]).astype(int)]
                # Should get an array of maxMOI**2 sums
                / np.sum(
                    frequencies_RR[1, y, : frequencies_RR.astype(int)[0, y, 0]]
                    * dvect[
                        correction_distance_matrix[  # TODO: Make sure multiplications are down the right axis (I believe they default to down columns, which is what I want)
                            y,
                            : frequencies_RR.astype(int)[
                                0, y, 0
                            ],  # TODO: Figure out how to do this more cleanly? (R impementation just used ":", assumed array had correct dimensions)
                            allrecrf[x, y, : maxMOI ** 2].astype(int),
                        ]
                    ],
                    axis=1,
                ),  # TODO: Verify it's the right axis?
            )
        likelihoodratio[x] = np.exp(np.sum(np.log(id_means)))

    z = np.random.uniform(size=nids)
    newclassification = classification
    newclassification[np.logical_and(classification == 0, z < likelihoodratio)] = 1
    newclassification[np.logical_and(classification == 1, z < 1 / likelihoodratio)] = 0
    classification = newclassification

    # propose new hidden states
    """
    # TODO: What does switch_hidden do? Is it entirely side effects? (Also,: can't run this yet, still waiting on implementation)
    for i in range(nids):
        switch_hidden(i)
    """

    # propose q (beta distribution is conjugate distribution for binomial process)
    q_prior_alpha = 0
    q_prior_beta = 0
    q_posterior_alpha = (
        q_prior_alpha + np.nansum(hidden0 == 1) + np.nansum(hiddenf == 1)
    )
    q_posterior_beta = q_prior_beta + np.nansum(hidden0 == 0) + np.nansum(hiddenf == 0)
    if q_posterior_alpha == 0:
        q_posterior_alpha = 1
    if q_posterior_beta == 0:  # TODO: Added this due to numpy warning, possibly remove?
        q_posterior_beta = 1
    qq = np.random.beta(q_posterior_alpha, q_posterior_beta)


@pytest.mark.xfail(reason='Not implemented in refactor')
def test_update_dvect():
    nids = 6
    nloci = 7

    classification = np.repeat(0, nids)
    classification[0] = 1
    mindistance = np.zeros((nids, nloci))

    dvect = np.ones(133)

    # Inputs needed: classification, mindistance, dvect
    if np.sum(classification == 1) >= 1:
        d_prior_alpha = 0
        d_prior_beta = 0
        d_posterior_alpha = d_prior_alpha + mindistance[classification == 1, :].size
        d_posterior_beta = d_prior_beta + np.sum(
            np.round(mindistance[classification == 1, :])
        )
        if d_posterior_beta == 0:
            d_posterior_beta = np.sum(mindistance[classification == 1, :])
        if (
            d_posterior_beta == 0
        ):  ## algorithm will get stuck if dposterior is allowed to go to 1 (TODO: Wait, so why is it setting d_posterior_beta to 1??)
            d_posterior_beta = 1

        dposterior = np.random.beta(d_posterior_alpha, d_posterior_beta)
        dvect = dposterior * (np.array(1 - dposterior) ** np.arange(0, dvect.size))
        dvect = dvect / np.sum(dvect)

    # TODO: Assert that dvect is updated correctly w/ beta function?


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
