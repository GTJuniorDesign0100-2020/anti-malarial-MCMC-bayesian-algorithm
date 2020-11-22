from typing import List

import numpy as np
import pandas as pd
import scipy.stats as sp_stats

from api.define_alleles import define_alleles
from api.findposteriorfrequencies import findposteriorfrequencies
import api.recrudescence_utils as recrudescence_utils
from api.site_instance_state import SiteInstanceState, SampleType, HiddenAlleleType
from api.switch_hidden import switch_hidden


class SavedState:
    def __init__(self, num_records, ids, classification, alleles0, allelesf, parameters):
        self.num_records = num_records
        self.ids = ids
        self.classification = classification
        self.alleles0 = alleles0
        self.allelesf = allelesf
        self.parameters = parameters
        self.posterior_df = None
        self.summary_stats_df = None

    def update_saved_state(
        self,
        state: SiteInstanceState,
        locinames: np.ndarray,
        iteration: int,
        burnin: int,
        record_interval: int):
        '''
        Updates the saved state instance
        '''
        if iteration <= burnin or iteration % record_interval != 0:
            return
        record_index = int((iteration - burnin) / record_interval)

        self.classification[:, record_index] = state.classification
        self.alleles0[:, :, record_index] = state.alleles0
        self.allelesf[:, :, record_index] = state.allelesf

        num_loci = locinames.size
        self.parameters[0, record_index] = state.qq
        self.parameters[1, record_index] = state.dposterior
        self.parameters[2 : (2 + num_loci), record_index] = state.frequencies_RR[1].max(axis=1)
        self.parameters[2 + num_loci : (2 + 2 * num_loci), record_index] = np.sum(
            state.frequencies_RR[1][:num_loci, :] ** 2
        )

    def post_process_saved_state(self, locinames: np.ndarray, num_ids: int, max_MOI: int):
        '''
        TODO: Verify correctness here
        Remove NaN columns from the saved state arrays, then calculate summary
        statistics for the algorithm run
        '''
        self.parameters = self.parameters[
            :, ~np.isnan(np.sum(self.parameters, axis=0))]
        self.classification = self.classification[
            :, ~np.isnan(np.sum(self.classification, axis=0))]

        self.posterior_df, self.summary_stats_df = self.get_summary_stats(locinames, num_ids, max_MOI)

    def get_summary_stats(self, locinames: np.ndarray, num_ids: int, max_MOI: int):
        '''
        Returns dataframes containing the summary statistics of the algorithm
        run
        :return: A tuple of (posterior_df, summary_stats_df)
        '''
        num_loci = locinames.size
        # TODO: What does this mean?
        # Find the mode of the hidden alleles on the first/last day
        modealleles = np.zeros((2 * num_ids, max_MOI * num_loci))
        for i in range(num_ids):
            for j in range(num_loci):
                modealleles[2 * i, j * max_MOI: (j + 1) * max_MOI] = sp_stats.mode(
                    self.alleles0[i, j * max_MOI: (j + 1) * max_MOI, :], axis=1
                )[0].ravel()

                modealleles[2 * i + 1, j * max_MOI: (j + 1) * max_MOI] = sp_stats.mode(
                    self.allelesf[i, j * max_MOI: (j + 1) * max_MOI, :], axis=1
                )[0].ravel()

        # TODO: Combined what? Why is this doubled?
        temp_combined = np.repeat(np.mean(self.classification, axis=1)[:num_ids], 2)
        # Reshape into a single column
        temp_combined = temp_combined.reshape(2 * num_ids, 1)
        posterior_matrix = np.concatenate((temp_combined, modealleles), axis=1)
        posterior_matrix_columns = [
            [f"{locus}_{i+1}" for i in range(max_MOI)] for locus in locinames
        ]
        posterior_matrix_columns = np.array(posterior_matrix_columns).flatten().tolist()
        posterior_matrix_columns.insert(0, "Prob Rec")
        posterior_df = pd.DataFrame(posterior_matrix, columns=posterior_matrix_columns)

        # TODO: Understand what's being printed here? Separate out into its own
        # function?
        summary_statisticsmatrix = np.concatenate((
                np.mean(self.parameters, axis=1).reshape(-1, 1),
                np.quantile(self.parameters, (0.25, 0.75), axis=1).T),
            axis=1)
        summary_statisticsmatrix = np.concatenate((
                summary_statisticsmatrix,
                np.append(
                    np.quantile(self.parameters[2 + num_loci:, :], (0.25, 0.75)),
                    np.mean(self.parameters[2 + num_loci:, :]),
                ).reshape(1, -1)))
        summary_statisticsmatrix = np.array([
                f"{summary_statisticsmatrix[i,0]:.2f} ({summary_statisticsmatrix[i,1]:.2f}, {summary_statisticsmatrix[i,2]:.2f})"
                for i in range(summary_statisticsmatrix.shape[0])
            ])
        summary_statistics_df = pd.DataFrame(
            summary_statisticsmatrix,
            index=["q", "d", *locinames.tolist(), *locinames.tolist(), "Mean diversity"])

        return posterior_df, summary_statistics_df


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
        seed: int=None) -> SavedState:
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
        :return: Returns the saved state (including the summary statistics) of
        the algorithm run
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
            self.ids,
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
            self.saved_state.update_saved_state(self.state, self.locinames,
                i, burnin, record_interval)
            print(f'MCMC Iteration {i + 1}')

        self.saved_state.post_process_saved_state(
            self.locinames,
            self.ids.size,
            self.max_MOI)

        return self.saved_state

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
        ids: np.ndarray,
        num_loci: int,
        max_MOI: int) -> SavedState:
        '''
        TODO: Reduce the number of parameters here
        Return a newly initialized SavedState object for a new algorithm run
        '''
        num_records = int((nruns - burnin) / record_interval)
        alleles_shape = (ids.size, max_MOI * num_loci, num_records)
        return SavedState(
            num_records=num_records,
            ids=ids,
            classification=recrudescence_utils.nan_array(
                (ids.size, num_records)),
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
        likelihood_ratios = cls._likelihood_ratios(
            state, num_ids, num_loci, max_MOI)

        cls._update_classifications(state, likelihood_ratios, num_ids, rand)

        # propose new hidden states
        # TODO: What does switch_hidden do? Is it entirely side effects?
        for i in range(num_ids):
            switch_hidden(i, num_loci, max_MOI, alleles_definitions_RR, state, rand)

        cls._update_q(state, rand)
        cls._update_dvect(state, rand)
        cls._update_frequencies(state, num_loci, max_MOI, rand)

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
            if id_means.all() != 0:
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
            # Should get an array of max_MOI**2 sums
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

    @classmethod
    def _update_classifications(
        cls,
        state: SiteInstanceState,
        likelihood_ratios: np.ndarray,
        num_ids: int,
        rand: np.random.RandomState):
        '''
        Update the recrudescence/reinfection classification of each sample,
        based on the current state's calculate likelihood ratios

        :param state: The current state variables of the algorithm
        :param likelihood_ratios:
        :param num_ids:
        :param rand: The random number generator to use
        :return: The new classifications (will also modify the state
        classifications)
        '''
        z = rand.uniform(size=num_ids)
        new_classifications = np.copy(state.classification)
        new_classifications[np.logical_and(
            state.classification == SampleType.REINFECTION.value,
            z < likelihood_ratios)] = SampleType.RECRUDESCENCE.value
        # Add slight offset so likelihood ratios don't give div by 0 error
        new_likelihood_ratios = likelihood_ratios + 1e-9
        new_classifications[np.logical_and(
            state.classification == SampleType.RECRUDESCENCE.value,
            z < 1.0 / new_likelihood_ratios)] = SampleType.REINFECTION.value
        state.classification = new_classifications
        return state.classification

    @classmethod
    def _update_q(cls, state: SiteInstanceState, rand: np.random.RandomState):
        '''
        TODO: Possibly move this to the state object itself?
        Propose a new value for q (the proportion of alleles that are hidden/
        not directly observed) and update it appropriately, based on the
        current state

        :param state: The current state variables of the algorithm
        :param rand: The random number generator to use
        :return: The updated q value (will also modify state.q)
        '''
        # TODO: What are the alpha/beta used for, in high-level terms? They seem
        # to be counts of observed/missing alleles in the hidden state?
        q_prior_alpha = 0
        q_prior_beta = 0

        q_posterior_alpha = (
            q_prior_alpha
            + np.nansum(state.hidden0 == HiddenAlleleType.MISSING.value)
            + np.nansum(state.hiddenf == HiddenAlleleType.MISSING.value))
        q_posterior_beta = (
            q_prior_beta
            + np.nansum(state.hidden0 == HiddenAlleleType.OBSERVED.value)
            + np.nansum(state.hiddenf == HiddenAlleleType.OBSERVED.value))

        # Edge case if there are no missing/observed alleles, to avoid div by 0
        if q_posterior_alpha == 0:
            q_posterior_alpha = 1
        if q_posterior_beta == 0:  # TODO: Added this due to numpy warning, possibly remove?
            q_posterior_beta = 1

        # propose new q (beta distribution is conjugate distribution for
        # binomial process)
        state.qq = rand.beta(q_posterior_alpha, q_posterior_beta)
        return state.qq

    @classmethod
    def _update_dvect(cls, state: SiteInstanceState, rand: np.random.RandomState):
        '''
        Update the distance vector (approximated using the geometric
        distribution) and dposterior based on the new state

        :param state: The current state variables of the algorithm
        :param rand: The random number generator to use
        :return: The modified dvector (will also modify state.dvect and
        state.dposterior)
        '''
        # only update if there is at least 1 recrudescent infection
        if np.sum(state.classification == SampleType.RECRUDESCENCE.value) == 0:
            return
        d_prior_alpha = 0
        d_prior_beta = 0
        min_recrudescence_distances = state.mindistance[state.classification == SampleType.RECRUDESCENCE.value, :]

        d_posterior_alpha = d_prior_alpha + min_recrudescence_distances.size
        d_posterior_beta = d_prior_beta + np.sum(
            np.round(min_recrudescence_distances))
        if d_posterior_beta == 0:
            d_posterior_beta = np.sum(min_recrudescence_distances)
        if d_posterior_beta == 0:  # algorithm will get stuck if dposterior is allowed to go to 1 (TODO: Wait, so why is it setting d_posterior_beta to 1??)
            d_posterior_beta = 1

        # TODO: Verify how this update actually works?
        state.dposterior = rand.beta(d_posterior_alpha, d_posterior_beta)
        #  update dvect (approximate using geometric distribution)
        state.dvect = state.dposterior * (
            np.array(1 - state.dposterior)**np.arange(0, state.dvect.size))
        state.dvect = state.dvect / np.sum(state.dvect)

        return state.dvect

    @classmethod
    def _update_frequencies(cls, state: SiteInstanceState, num_loci: int, max_MOI: int, rand: np.random.RandomState):
        '''
        # TODO: Verify what this does, in full?
        Update the (allele?) frequencies based on the current algorithm state
        '''
        # first, remove recrudescing alleles from calculations
        tempdata = np.copy(state.recoded0)
        # TODO: Ask about line: "sapply(which(classification == 1), function
        # (x) tempdata[x,recr0[x,]] <<- 0)" in the original code.
        recrudescent_allele_indices = np.where(state.classification == SampleType.RECRUDESCENCE.value)[0]  # [0] unboxes the tuple.

        for recrud_index in recrudescent_allele_indices:
            columns = state.recr0[recrud_index, :].astype(int)
            tempdata[recrud_index, columns] = 0

        # actually update the frequencies
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
                state.frequencies_RR,
                rand)
