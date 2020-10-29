from collections import namedtuple
from typing import List

import numpy as np
import pandas as pd

from define_alleles import define_alleles
from findposteriorfrequencies import findposteriorfrequencies
import recrudescence_utils
from site_instance_state import SiteInstanceState
from switch_hidden import switch_hidden


SavedState = namedtuple(
    'SavedState',
    'num_records classification alleles0 allelesf parameters')


class AlgorithmSiteInstance:
    '''
    Handles running the MCMC malaria recrudescence algorithm for a single
    "arm"/site's data
    '''

    def __init__(
        self,
        genotypedata_RR: pd.DataFrame,
        additional_neutral: pd.DataFrame,
        locirepeats: List[int]):
        '''
        Sets up the initial data structures needed before running the algorithm

        :param genotypedata_RR: The dataset with day 0/day of failure data for
        each sample in the site
        :param additional_neutral: The set of additional initial samples taken
        to provide background information about allele prevalence
        '''
        # MOI = multiplicity of infection
        self.max_MOI = self._get_max_MOI(genotypedata_RR)

        # Get the unique Sample IDs/loci in the dataset
        # NOTE: pd.unique used instead of np.unique to preserve ordering
        self.ids = recrudescence_utils.get_sample_ids(genotypedata_RR, 'Day 0')
        self.locinames = pd.unique(genotypedata_RR.columns[1:].str.split("_").str[0])

        # TODO: Should this be here or on the state?
        self.alleles_definitions_RR = self._get_allele_definitions(
            genotypedata_RR, additional_neutral, self.locinames.size, locirepeats)

        # Set up the initial state for the algorithm
        self.state = SiteInstanceState(
            self.ids,
            self.locinames,
            self.max_MOI,
            genotypedata_RR,
            additional_neutral,
            self.alleles_definitions_RR
        )

    def run_algorithm(
        self,
        jobname: str,
        nruns: int=1000,
        burnin: int=100,
        record_interval: int=10,
        seed: int=None):
        '''
        Runs the actual algorithm on this site's data and returns the results
        TODO: Expand on this

        :param jobname: The name to label output data/files from this run with
        :param nruns: (optional) The number of iterations to run the algorithm;
        more iterations will take longer, but will be more accurate
        :param burnin: (optional) The number of initial iterations to not store
        data for (as the algorithm won't have converged/have usable results by
        that point)
        :param record_interval: (optional) Record data from the algorithm every
        "record_interval" iterations
        :param seed: (optional) The seed to use for random numbers when running
        (defaults to completely random)
        :return: TODO: Currently returns the saved classifications, parameter
        matrix, and ids from this run (determine if more/less should be returned?)
        '''
        self.state.randomize_initial_assignments(
            self.ids.size,
            self.locinames.size,
            self.max_MOI,
            self.alleles_definitions_RR,
            seed)

        self.saved_state = self._get_initial_saved_state(
            nruns,
            burnin,
            record_interval,
            self.ids.size,
            self.locinames.size,
            self.max_MOI)

        rand = np.random.RandomState(seed=seed)
        for i in range(nruns):
            self._run_mcmc(
                self.state,
                self.alleles_definitions_RR,
                self.ids.size,
                self.locinames.size,
                self.max_MOI,
                rand)
            self._update_saved_state(self.state, self.saved_state,
                i, burnin, record_interval)
            print(f'MCMC Iteration {i + 1}')

        # TODO: Add in summary stats

        return (
            self.saved_state.classification,
            self.saved_state.parameters,
            self.ids
        )

    @classmethod
    def _get_max_MOI(cls, genotypedata_RR: pd.DataFrame) -> int:
        '''
        Get the maximum "multiplicity of infection" (MOI) in the dataset (i.e.
        the maximum number of strains there are for a single locus in the
        dataset)
        '''
        return int(np.nanmax(
            pd.to_numeric(genotypedata_RR.columns.str.split('_').str[1])
        ))

    @classmethod
    def _get_allele_definitions(
        cls,
        genotypedata_RR: pd.DataFrame,
        additional_neutral: pd.DataFrame,
        num_loci: int,
        locirepeats: List[int]) -> pd.DataFrame:
        '''
        Get the allele definitions in the dataset (i.e. the set of ranges/
        "bins" the allele lengths can fall into)
        '''
        maxalleles = 30
        k = np.repeat(maxalleles, num_loci)
        alleles_definitions = define_alleles(
            pd.concat([genotypedata_RR, additional_neutral]), locirepeats, k
        )
        return alleles_definitions

    @classmethod
    def _get_initial_saved_state(
        cls,
        nruns: int,
        burnin: int,
        record_interval: int,
        num_ids: int,
        num_loci: int,
        max_MOI: int) -> SavedState:
        '''
        TODO: Reduce the number of parameters here
        Return a newly initialized SavedState object for a new algorithm run
        '''
        num_records = int((nruns - burnin) / record_interval)
        alleles_shape = (num_ids, max_MOI * num_loci, num_records)
        return SavedState(
            num_records=num_records,
            classification=recrudescence_utils.nan_array(
                (num_ids, num_records)),
            alleles0=recrudescence_utils.nan_array(alleles_shape),
            allelesf=recrudescence_utils.nan_array(alleles_shape),
            parameters=recrudescence_utils.nan_array(
                (2 + 2 * num_loci, num_records))
        )

    @classmethod
    def _run_mcmc(
        cls,
        state: SiteInstanceState,
        alleles_definitions_RR: List[pd.DataFrame],
        num_ids: int,
        num_loci: int,
        max_MOI: int,
        rand: np.random.RandomState):
        '''
        NOTE: Modifies state via side-effects for performance reasons (no
        return)
        TODO: Cut down on parameters
        TODO: Explain this more fully
        Runs 1 iteration of the main MCMC algorithm loop, updating the state
        appropriately
        '''
        # propose new classification
        likelihoodratio = cls._likelihood_ratios(
            state, num_ids, num_loci, max_MOI)

        z = rand.uniform(size=num_ids)
        newclassification = state.classification
        newclassification[np.logical_and(state.classification == 0, z < likelihoodratio)] = 1
        newclassification[np.logical_and(state.classification == 1, z < 1 / likelihoodratio)] = 0
        state.classification = newclassification

        # propose new hidden states
        # TODO: What does switch_hidden do? Is it entirely side effects?
        for i in range(num_ids):
            switch_hidden(i, num_loci, max_MOI, alleles_definitions_RR, state)

        # propose q (beta distribution is conjugate distribution for binomial process)
        q_prior_alpha = 0
        q_prior_beta = 0
        q_posterior_alpha = (
            q_prior_alpha + np.nansum(state.hidden0 == 1) + np.nansum(state.hiddenf == 1)
        )
        q_posterior_beta = q_prior_beta + np.nansum(state.hidden0 == 0) + np.nansum(state.hiddenf == 0)
        if q_posterior_alpha == 0:
            q_posterior_alpha = 1
        if q_posterior_beta == 0:  # TODO: Added this due to numpy warning, possibly remove?
            q_posterior_beta = 1
        state.qq = rand.beta(q_posterior_alpha, q_posterior_beta)

        #  update dvect (approximate using geometric distribution)
        # only if there is at least 1 recrudescent infection
        if np.sum(state.classification == 1) >= 1:
            d_prior_alpha = 0
            d_prior_beta = 0
            d_posterior_alpha = d_prior_alpha + state.mindistance[state.classification == 1, :].size
            d_posterior_beta = d_prior_beta + np.sum(
                np.round(state.mindistance[state.classification == 1, :])
            )
            if d_posterior_beta == 0:
                d_posterior_beta = np.sum(state.mindistance[state.classification == 1, :])
            if (
                d_posterior_beta == 0
            ):  ## algorithm will get stuck if dposterior is allowed to go to 1 (TODO: Wait, so why is it setting d_posterior_beta to 1??)
                d_posterior_beta = 1

            state.dposterior = np.random.beta(d_posterior_alpha, d_posterior_beta)
            state.dvect = state.dposterior * (
                np.array(1 - state.dposterior)**np.arange(0, state.dvect.size))
            state.dvect = state.dvect / np.sum(state.dvect)

        # update frequencies
        # # first, remove recrudescing alleles from calculations
        # TODO: Deep or Shallow copy?
        tempdata = np.copy(state.recoded0)
        # TODO: Ask about this line *** sapply(which(classification == 1), function (x) tempdata[x,recr0[x,]] <<- 0) *** in the OG code.
        recrudescent_alleles = np.where(state.classification == 1)[0] # unboxes the tuple.
        end = state.recr0[recrudescent_alleles, :]
        if (len(recrudescent_alleles) > 0):
            for recrud in recrudescent_alleles:
                columns = state.recr0[recrud,].astype(int)
                tempdata[recrud, columns] = 0

        tempdata = np.concatenate((tempdata, state.recodedf), axis=0)
        for i in range(num_loci):
            # TODO: Verify ignoring recoded_additional_neutral is correct if it
            # wasn't initialized?
            tempdata_with_recoded = tempdata
            if state.recoded_additional_neutral.size > 0:
                tempdata_with_recoded = np.concatenate((tempdata, state.recoded_additional_neutral), axis=0)
            findposteriorfrequencies(
                i,
                tempdata_with_recoded,
                max_MOI,
                state.frequencies_RR)

    @classmethod
    def _likelihood_ratios(
        cls,
        state: SiteInstanceState,
        num_ids: int,
        num_loci: int,
        max_MOI: int):
        '''
        TODO: Explain what this does?
        Returns the likelihood ratio for each sample in the dataset, given our
        current state
        '''
        likelihoodratio = np.zeros(num_ids)
        # TODO: Finish vectorizing this

        for x in range(num_ids):
            # id mean for what?
            id_means = np.zeros(num_loci)
            for y in range(num_loci):
                id_means[y] = np.nanmean(
                    cls._likelihood_inner_loop(state, max_MOI, x, y)
                )
            likelihoodratio[x] = np.exp(np.sum(np.log(id_means)))
        return likelihoodratio

    @classmethod
    def _likelihood_inner_loop(cls, state, max_MOI, x: int, y: int):
        '''
        TODO: What does this actually do?
        Returns a 1D vector of max_MOI**2 length
        '''
        def non_nan(array: np.ndarray):
            # Needed since converting NaN to int has undefined behavior
            return array[~np.isnan(array)]

        dvect_indices = np.round(non_nan(state.alldistance[x, y, :])).astype(int)
        return (state.dvect[dvect_indices] /
            # Should get an array of maxMOI**2 sums
            np.sum(
                # TODO: Make sure multiplications are down the right axis (I believe each element in the frequencies_RR 1D vector should multiply across 1 dvect row)
                # Double-transpose to multiply across rows, not columns
                (state.frequencies_RR[1][y, :int(state.frequencies_RR[0][y])]
                * state.dvect[
                    state.correction_distance_matrix[y][
                        :,
                        non_nan(state.allrecrf[x, y, :max_MOI**2]).astype(int),
                    ].astype(int)
                ].T).T,
                axis=0
            ))

    def _update_saved_state(
        self,
        state: SiteInstanceState,
        saved_state: SavedState,
        iteration: int,
        burnin: int,
        record_interval: int):
        '''
        Updates the saved state instance
        '''
        if iteration <= burnin or iteration % record_interval != 0:
            return
        record_index = int((iteration - burnin) / record_interval)

        saved_state.classification[:, record_index] = state.classification
        saved_state.alleles0[:, :, record_index] = state.alleles0
        saved_state.allelesf[:, :, record_index] = state.allelesf

        num_loci = self.locinames.size
        saved_state.parameters[0, record_index] = state.qq
        saved_state.parameters[1, record_index] = state.dposterior
        saved_state.parameters[2 : (2 + num_loci), record_index] = state.frequencies_RR[1].max(axis=1)
        saved_state.parameters[2 + num_loci : (2 + 2 * num_loci), record_index] = np.sum(
            state.frequencies_RR[1][:num_loci, :] ** 2
        )

    @classmethod
    def _remove_nan_from_array(cls, array: np.ndarray) -> np.ndarray:
        '''
        Returns a version of the array with NaNs removed (does not modify the
        original array)

        :param array: The array to process
        :return: A version of the array with all np.nan values removed
        '''
        return array[:, np.sum(~np.isnan(array), axis=1)]
