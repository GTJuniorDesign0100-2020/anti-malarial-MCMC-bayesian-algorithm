from collections import namedtuple
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
from api.recrudescence_file_parser import RecrudescenceFileParser
from api.switch_hidden import switch_hidden


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

# TODO: Copy-pasted from state setup in test_algorithm_site_instance (w/ minor tweaks)
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

    state.dvect = np.zeros(300)
    state.dvect[:3] = np.array([0.75, 0.2, 0.05])
    state.MOI0 = np.ones(num_ids, dtype=int)
    state.MOIf = np.ones(num_ids, dtype=int)
    state.hidden0 = state.hidden0[:num_ids, :max_MOI*num_loci]
    state.hiddenf = state.hiddenf[:num_ids, :max_MOI*num_loci]
    state.alleles0 = np.array([
        [223.4, 262, 232, 0, 0, 103.7],
        [229.3, 250, 0, 0, 0, 139.1],
        [221.3, 250, 222, 0, 0, 151.6],
    ])
    state.allelesf = np.array([
        [225.5, 222, 0, 0, 0, 124.1],
        [243.6, 224, 0, 0, 0, 139.0],
        [231.5, 250, 0, 0, 0, 145.0],
    ])
    state.recoded0 = np.array([
        [5, 7, 8, 0, 0, 3],
        [12, 7, 0, 0, 0, 2],
        [2, 10, 9, 0, 0, 13],
    ])
    state.recodedf = np.array([
        [6, 1, 0, 0, 0, 1],
        [3, 5, 0, 0, 0, 10],
        [7, 2, 0, 0, 0, 6],
    ])
    state.recr0 = np.array([
        [0, 2, 4],
        [1, 3, 5],
        [0, 2, 4],
    ])
    state.recrf = np.array([
        [1, 3, 5],
        [0, 2, 4],
        [0, 2, 4],
    ])

    return BasicState(state, num_ids, num_loci, max_MOI)

def test_switch_hidden_allele_reinfection(mcmc_initial_state):
    state = mcmc_initial_state.state
    num_ids = mcmc_initial_state.num_ids
    num_loci = mcmc_initial_state.num_loci
    max_MOI = mcmc_initial_state.max_MOI

    # Copy parameters switch_hidden modifies so the inputs don't change
    copy_alleles0 = state.alleles0.copy()
    copy_allelesf = state.allelesf.copy()
    copy_recoded0 = state.recoded0.copy()
    copy_recodedf = state.recodedf.copy()
    copy_mindistance = state.mindistance.copy()
    copy_alldistance = state.alldistance.copy()
    copy_recr0 = state.recr0.copy()
    copy_recrf = state.recrf.copy()
    copy_allrecrf = state.allrecrf.copy()

    accumulate_alleles0 = np.zeros((num_ids, max_MOI*num_loci))
    accumulate_allelesf = np.zeros((num_ids, max_MOI*num_loci))
    accumulate_mindist = np.zeros((num_ids, num_loci))

    rand = np.random.RandomState(2020)
    num_iterations = 1000
    for i in range(num_iterations):
        state.alleles0 = copy_alleles0.copy()
        state.allelesf = copy_allelesf.copy()
        state.recoded0 = copy_recoded0.copy()
        state.recodedf = copy_recodedf.copy()
        state.mindistance = copy_mindistance.copy()
        state.alldistance = copy_alldistance.copy()
        state.recr0 = copy_recr0.copy()
        state.recrf = copy_recrf.copy()
        state.allrecrf = copy_allrecrf.copy()

        switch_hidden(
            x=0,
            nloci=num_loci,
            maxMOI=max_MOI,
            alleles_definitions_RR=alleles_definitions_RR,
            state=state,
            rand=rand)

        accumulate_alleles0 += state.alleles0
        accumulate_allelesf += state.allelesf
        accumulate_mindist += state.mindistance

    accumulate_alleles0 /= num_iterations
    accumulate_allelesf /= num_iterations
    accumulate_mindist /= num_iterations

    # Verify only the 1st element has updated at this point
    np.testing.assert_almost_equal(accumulate_alleles0[0, 1:], copy_alleles0[0, 1:])
    np.testing.assert_almost_equal(accumulate_alleles0[1:], copy_alleles0[1:])
    assert 0.1 <= accumulate_alleles0[0, 0] <= 15.0
    np.testing.assert_almost_equal(accumulate_allelesf, copy_allelesf)
    np.testing.assert_almost_equal(accumulate_mindist, copy_mindistance)


def test_switch_hidden_allele_recrudescence(mcmc_initial_state):
    state = mcmc_initial_state.state
    num_ids = mcmc_initial_state.num_ids
    num_loci = mcmc_initial_state.num_loci
    max_MOI = mcmc_initial_state.max_MOI

    # Copy parameters switch_hidden modifies so the inputs don't change
    copy_alleles0 = state.alleles0.copy()
    copy_allelesf = state.allelesf.copy()
    copy_recoded0 = state.recoded0.copy()
    copy_recodedf = state.recodedf.copy()
    copy_mindistance = state.mindistance.copy()
    copy_alldistance = state.alldistance.copy()
    copy_recr0 = state.recr0.copy()
    copy_recrf = state.recrf.copy()
    copy_allrecrf = state.allrecrf.copy()

    accumulate_alleles0 = np.zeros((num_ids, max_MOI*num_loci))
    accumulate_allelesf = np.zeros((num_ids, max_MOI*num_loci))
    accumulate_mindist = np.zeros((num_ids, num_loci))

    rand = np.random.RandomState(1337)
    num_iterations = 1000
    for i in range(num_iterations):
        state.alleles0 = copy_alleles0.copy()
        state.allelesf = copy_allelesf.copy()
        state.recoded0 = copy_recoded0.copy()
        state.recodedf = copy_recodedf.copy()
        state.mindistance = copy_mindistance.copy()
        state.alldistance = copy_alldistance.copy()
        state.recr0 = copy_recr0.copy()
        state.recrf = copy_recrf.copy()
        state.allrecrf = copy_allrecrf.copy()

        switch_hidden(
            x=2,
            nloci=num_loci,
            maxMOI=max_MOI,
            alleles_definitions_RR=alleles_definitions_RR,
            state=state,
            rand=rand)

        accumulate_alleles0 += state.alleles0
        accumulate_allelesf += state.allelesf
        accumulate_mindist += state.mindistance

    accumulate_alleles0 /= num_iterations
    accumulate_allelesf /= num_iterations
    accumulate_mindist /= num_iterations

    # Verify only the 1st element has updated at this point
    np.testing.assert_almost_equal(accumulate_alleles0[2, 1:], copy_alleles0[2, 1:])
    np.testing.assert_almost_equal(accumulate_alleles0[:2], copy_alleles0[:2])
    assert 0.1 <= accumulate_alleles0[2, 0] <= 15.0
    np.testing.assert_almost_equal(accumulate_allelesf, copy_allelesf)
    np.testing.assert_almost_equal(accumulate_mindist, copy_mindistance)
