from collections import namedtuple
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

    return BasicState(state, num_ids, num_loci, max_MOI)


def test_likelihood_ratio_inner_loop_initial(mcmc_initial_state):
    expected_inner_values = np.array([8.265306, 8.265306, 8.265306, 8.265306])

    inner_values = AlgorithmSiteInstance._likelihood_inner_loop(
        mcmc_initial_state.state,
        mcmc_initial_state.max_MOI,
        0,
        0)

    np.testing.assert_array_almost_equal(inner_values, expected_inner_values, decimal=5)


def test_likelihood_ratio_inner_loop_middle(mcmc_initial_state):
    expected_inner_values = np.array([0.7777778, 0.7777778, 0.7777778, 0.7777778])

    inner_values = AlgorithmSiteInstance._likelihood_inner_loop(
        mcmc_initial_state.state,
        mcmc_initial_state.max_MOI,
        1,
        2)

    np.testing.assert_array_almost_equal(inner_values, expected_inner_values, decimal=5)


def test_likelihood_ratio(mcmc_initial_state):
    # NOTE: Middle number changes by +/-3 in the R code?
    expected_likelihood_ratios = np.array([56.024203, 1.438889, 0.0])

    likelihood_ratios = AlgorithmSiteInstance._likelihood_ratios(
        mcmc_initial_state.state,
        mcmc_initial_state.num_ids,
        mcmc_initial_state.num_loci,
        mcmc_initial_state.max_MOI)

    np.testing.assert_array_almost_equal(likelihood_ratios, expected_likelihood_ratios, decimal=5)
