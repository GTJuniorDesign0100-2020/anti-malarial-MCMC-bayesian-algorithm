from collections import namedtuple
import cProfile
import os
import sys

import numpy as np
import pandas as pd
import pytest
import scipy.stats as sp_stats

# Add parent directory to search path, so we can import those files
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.algorithm_instance import AlgorithmInstance
from api.algorithm_site_instance import AlgorithmSiteInstance, SiteInstanceState
from api.findposteriorfrequencies import findposteriorfrequencies
from api.recrudescence_file_parser import RecrudescenceFileParser
import api.recrudescence_utils as recrudescence_utils


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


BasicState = namedtuple('BasicState', 'state num_ids num_loci max_MOI')

@pytest.fixture
def mcmc_initial_state():
    state = SiteInstanceState(expected_ids, expected_locinames, expected_maxMOI, genotypedata_RR, additional_neutral, alleles_definitions_RR)
    state.randomize_initial_assignments(
        expected_ids.size,
        expected_locinames.size,
        expected_maxMOI,
        alleles_definitions_RR
    )
    # NOTE: Replaces the randomized arrays with hardcoded data. This is entirely
    # stubbed data, and may not reflect actual data at all
    num_ids = 3
    num_loci = 3
    max_MOI = 2
    state.alldistance = np.arange(0, 3.6, step=0.1).reshape(3, 3, max_MOI**2)
    state.mindistance = np.array([
        [0.0, 0.4, 0.8],
        [1.2, 1.6, 2.0],
        [2.4, 2.8, 3.2],
    ])
    state.allrecrf = np.floor(np.arange(0, 9, step=0.25).reshape(3, 3, max_MOI**2))
    state.classification = np.array([0, 1, 1])

    return BasicState(state, num_ids, num_loci, max_MOI)


def test_likelihood_ratio_inner_loop_initial(mcmc_initial_state):
    expected_inner_values = np.array([8.265306, 8.265306, 8.265306, 8.265306])

    inner_values = AlgorithmSiteInstance._likelihood_inner_loop(
        mcmc_initial_state.state,
        0)

    np.testing.assert_array_almost_equal(inner_values[0], expected_inner_values, decimal=5)


def test_likelihood_ratio_inner_loop_middle(mcmc_initial_state):
    expected_inner_values = np.array([0.7777778, 0.7777778, 0.7777778, 0.7777778])

    inner_values = AlgorithmSiteInstance._likelihood_inner_loop(
        mcmc_initial_state.state,
        2)

    np.testing.assert_array_almost_equal(inner_values[1], expected_inner_values, decimal=5)


def test_likelihood_ratio(mcmc_initial_state):
    expected_likelihood_ratios = np.array([56.024203, 1.438889, 0.0])

    likelihood_ratios = AlgorithmSiteInstance._likelihood_ratios(
        mcmc_initial_state.state,
        mcmc_initial_state.num_ids,
        mcmc_initial_state.num_loci)

    np.testing.assert_array_almost_equal(likelihood_ratios, expected_likelihood_ratios, decimal=5)


def test_likelihood_ratio_with_nans(mcmc_initial_state):
    expected_likelihood_ratios = np.array([88.459268, 1.438889, 0.0])

    mcmc_initial_state.state.alldistance[0, 1, 2:] = np.nan
    mcmc_initial_state.state.allrecrf[0, 1, 2:] = np.nan
    likelihood_ratios = AlgorithmSiteInstance._likelihood_ratios(
        mcmc_initial_state.state,
        mcmc_initial_state.num_ids,
        mcmc_initial_state.num_loci)

    np.testing.assert_array_almost_equal(likelihood_ratios, expected_likelihood_ratios, decimal=5)


def test_likelihood_ratio_speed(mcmc_initial_state):
    cProfile.runctx('''for i in range(1000):
        AlgorithmSiteInstance._likelihood_ratios(
        mcmc_initial_state.state,
        mcmc_initial_state.num_ids,
        mcmc_initial_state.num_loci)''', globals(), locals(), 'new_likelihood.cprof')


def test_updating_classifications(mcmc_initial_state):
    expected_avg_classifications = np.array([1.0, 0.31, 0.0])

    original_classifications = np.copy(mcmc_initial_state.state.classification)
    classifications_sum = np.zeros(original_classifications.size)
    likelihood_ratios = np.array([56.024203, 1.438889, 0.0])
    rand = np.random.RandomState(2020)
    for i in range(1000):
        classifications = AlgorithmSiteInstance._update_classifications(
            mcmc_initial_state.state,
            likelihood_ratios,
            mcmc_initial_state.num_ids,
            rand)
        classifications_sum += classifications
        mcmc_initial_state.state.classification = original_classifications

    avg_classes = classifications_sum / 1000.0
    # TODO: Find a more robust way of testing this (since there's still a tiny chance it randomly falls outside this range?)
    np.testing.assert_array_almost_equal(avg_classes, expected_avg_classifications, decimal=2)


def test_updating_q(mcmc_initial_state):
    rand = np.random.RandomState(2020)
    q_sum = 0.0
    for i in range(1000):
        AlgorithmSiteInstance._update_q(mcmc_initial_state.state, rand)
        q_sum += mcmc_initial_state.state.qq

    q_average = q_sum / 1000.0
    # TODO: Find a more robust way of testing this (since there's still a tiny chance it randomly falls outside this range?)
    assert 0.375 <= q_average <= 0.380


def test_updating_dvect(mcmc_initial_state):
    rand = np.random.RandomState(2020)
    dposterior_sum = 0
    dvect_sum = np.zeros(133)
    for i in range(1000):
        AlgorithmSiteInstance._update_dvect(mcmc_initial_state.state, rand)
        dposterior_sum += mcmc_initial_state.state.dposterior
        dvect_sum += mcmc_initial_state.state.dvect

    dpost_average = dposterior_sum / 1000.0
    dvect_average = dvect_sum / 1000.0
    assert 0.31 <= dpost_average <= 0.32
    np.testing.assert_approx_equal(dvect_average[0], dpost_average)
    assert 0.202 <= dvect_average[1] <= 0.206


@pytest.mark.xfail(reason='Still need to figure out correct test output')
def test_updating_frequencies(mcmc_initial_state):
    rand = np.random.RandomState(2020)

    AlgorithmSiteInstance._update_frequencies(
        mcmc_initial_state.state,
        mcmc_initial_state.num_loci,
        mcmc_initial_state.max_MOI,
        rand)


@pytest.mark.xfail(reason='Still need to figure out correct test output')
def test_find_posterior_frequencies():
    # TODO: Frequencies don't update at all in the R code with this test case?
    # It clearly does over enough iterations, but...?
    expected_frequencies = [
        np.array([2, 4, 5]),
        np.array([
            [0.5, 0.5, 0, 0, 0],
            [0.25, 0.25, 0.25, 0.25, 0],
            [0.5, 0.125, 0.125, 0.125, 0.125]
        ]),
        np.array([0.0, 0.0, 0.167705])
    ]

    rand = np.random.RandomState(2020)
    test_frequencies = [
        np.array([2, 4, 5]),
        np.array([
            [0.5, 0.5, 0, 0, 0],
            [0.25, 0.25, 0.25, 0.25, 0],
            [0.5, 0.125, 0.125, 0.125, 0.125]
        ]),
        np.array([0.0, 0.0, 0.167705])
    ]
    test_tempdata = np.array([
        [0, 1, 0, 3, 0, 5],
        [1, 0, 2, 0, 4, 0],
        [3, 0, 0, 0, 0, 0],
    ])

    findposteriorfrequencies(0, test_tempdata, 2, test_frequencies, rand)

    np.testing.assert_array_equal(test_frequencies[0], expected_frequencies[0])
    np.testing.assert_array_almost_equal(test_frequencies[1], expected_frequencies[1])
    np.testing.assert_array_almost_equal(test_frequencies[2], expected_frequencies[2])


@pytest.mark.xfail(reason='Still need to figure out correct test output')
def test_updating_saved_state(mcmc_initial_state):
    rand = np.random.RandomState(2020)


@pytest.mark.xfail(reason='Still need to figure out correct test output')
def test_finding_allele_modes(mcmc_initial_state):
    rand = np.random.RandomState(2020)


@pytest.mark.xfail(reason='Still need to figure out correct test output')
def test_getting_summary_stats(mcmc_initial_state):
    rand = np.random.RandomState(2020)
