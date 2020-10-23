import enum
import typing

import numpy as np
import pandas as pd

from calculate_frequencies import calculate_frequencies3
from define_alleles import define_alleles
import mcmc
import recrudescence_utils


class SampleClass(enum.Enum):
    REINFECTION = 0
    RECRUDESCENCE = 1


class AlgorithmSiteInstance:
    '''
    Handles running the MCMC malaria recrudescence algorithm for a single
    "arm"/site's data
    '''

    def __init__(
        self,
        genotypedata_RR: pd.DataFrame,
        additional_neutral: pd.DataFrame,
        locirepeats: typing.List[int]):
        '''
        Sets up the initial data structures needed before running the algorithm
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



        # =====================================================================
        # TODO: Remove this!
        self.genotypedata_RR = genotypedata_RR
        self.additional_neutral = additional_neutral
        self.locirepeats = locirepeats
        # =====================================================================

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
        :return: TODO:
        '''
        self.state.randomize_initial_assignments(
            self.ids.size,
            self.locinames.size,
            self.max_MOI,
            self.alleles_definitions_RR,
            seed)
        # TODO: Actually implement this!
        return mcmc.onload(
            self.genotypedata_RR,
            self.additional_neutral,
            self.locirepeats,
            nruns,
            burnin,
            record_interval,
            jobname)

    @classmethod
    def _get_max_MOI(cls, genotypedata_RR: pd.DataFrame) -> int:
        '''
        Get the maximum "multiplicity of infection" (MOI) in the dataset
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
        locirepeats: typing.List[int]) -> pd.DataFrame:
        '''
        TODO: Elaborate
        Get the allele definitions in the dataset
        '''
        maxalleles = 30
        k = np.repeat(maxalleles, num_loci)
        alleles_definitions = define_alleles(
            pd.concat([genotypedata_RR, additional_neutral]), locirepeats, k
        )
        return alleles_definitions


class SiteInstanceState:
    '''
    A class that holds the current state for a single "arm"/site instance when
    running the malaria recrudescence algorithm
    '''

    def __init__(
        self,
        ids: np.ndarray,
        locinames: np.ndarray,
        maxMOI: int,
        genotypedata_RR: pd.DataFrame,
        additional_neutral: pd.DataFrame,
        alleles_definitions_RR: pd.DataFrame):
        '''
        Sets up the deterministic initial state of the algorithm
        TODO: Elaborate, shorten argument list?
        '''
        self._create_empty_state(ids.size, locinames.size, maxMOI)

        self.MOI0, self.MOIf = self._calculate_sample_MOI(
            genotypedata_RR, ids, locinames)
        self._initialize_alleles(
            genotypedata_RR, alleles_definitions_RR, locinames, maxMOI)

        # estimate frequencies (TODO: Unsure if this should be here, since it's
        # static state?)
        self.frequencies_RR = calculate_frequencies3(
            pd.concat([genotypedata_RR, additional_neutral]),
            alleles_definitions_RR)

    def _create_empty_state(self, num_ids: int, num_loci: int, max_MOI: int):
        '''
        Creates the initial empty data structures in the state to be populated
        later
        TODO: Explain what each created structure is and what it's used for
        '''
        self.alleles0 = np.zeros((num_ids, max_MOI * num_loci))
        self.recoded0 = np.zeros((num_ids, max_MOI * num_loci))
        self.hidden0 = np.full_like(np.empty((num_ids, max_MOI * num_loci)), np.nan)
        self.recr0 = np.full_like(np.empty((num_ids, num_loci)), np.nan)
        self.recr_repeats0 = np.full_like(
            np.empty((num_ids, num_loci)), np.nan
        )  # number of times recrudescing allele is repeated on day 0
        self.recr_repeatsf = np.full_like(
            np.empty((num_ids, num_loci)), np.nan
        )  # number of times recrudescing allele is repeated on day 0
        self.allelesf = np.zeros((num_ids, max_MOI * num_loci))
        self.recodedf = np.zeros((num_ids, max_MOI * num_loci))
        self.hiddenf = np.full_like(np.empty((num_ids, max_MOI * num_loci)), np.nan)
        self.recrf = np.full_like(np.empty((num_ids, num_loci)), np.nan)
        self.mindistance = np.zeros((num_ids, num_loci))
        self.alldistance = np.full_like(np.empty((num_ids, num_loci, max_MOI ** 2)), np.nan)
        self.allrecrf = np.full_like(np.empty((num_ids, num_loci, max_MOI ** 2)), np.nan)
        self.classification = np.repeat(SampleClass.REINFECTION, num_ids)

    @classmethod
    def _calculate_sample_MOI(
        cls,
        genotypedata_RR: pd.DataFrame,
        ids: np.ndarray,
        locinames: np.ndarray):
        '''
        TODO: Verify the details of what this is doing
        Calculates the "multiplicity of infections" (MOIs) on day 0 and the day
        of failure for each sample in the dataset. Returns these MOIs via 2
        arrays
        # TODO: Explain what multiplicity of infection actually means

        :param genotypedata_RR: The dataset holding the samples
        :param ids: The sample ids in the dataset to calculate MOIs for
        :locinames: The names of the unique loci in the dataset
        :return: A tuple of 1D numpy arrays of length "ids.size", (MOI0, MOIf)
        (which contain the multiplicities on the first/last day for each sample,
        respectively)
        '''
        MOI0 = np.repeat(0, ids.size)
        MOIf = np.repeat(0, ids.size)
        for i, ID in enumerate(ids):
            for lociname in locinames:
                locicolumns = genotypedata_RR.columns.str.contains(f"{lociname}_")

                num_alleles0 = np.count_nonzero(
                    ~genotypedata_RR.loc[
                        genotypedata_RR["Sample ID"].str.contains(f"{ID}_ Day 0"), locicolumns
                    ].isna()
                )
                num_allelesf = np.count_nonzero(
                    ~genotypedata_RR.loc[
                        genotypedata_RR["Sample ID"].str.contains(f"{ID}_ Day Failure"),
                        locicolumns,
                    ].isna()
                )

                MOI0[i] = np.max([MOI0[i], num_alleles0])
                MOIf[i] = np.max([MOIf[i], num_allelesf])
        return MOI0, MOIf

    def _initialize_alleles(
        self,
        genotypedata_RR: pd.DataFrame,
        alleles_definitions_RR: pd.DataFrame,
        locinames: np.ndarray,
        max_MOI: int):
        '''
        TODO: Verify what this actually does?
        Initialize the alleles/recoded state, filling them with the appropriate
        initial/failure data from the dataframe

        :param genotypedata_RR: The genotype data for this site the sample data
        :param alleles_definitions_RR: TODO:
        :param locinames: The names of the loci in genotypedata_RR to get
        alleles for
        :param max_MOI: The maximum multiplicity of infection in the genotype
        data
        '''
        for i, locus in enumerate(locinames):
            oldalleles, newalleles = self._get_original_alleles(
                genotypedata_RR, alleles_definitions_RR, locus, i)

            startColumn = max_MOI * i  # NOTE: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
            endColumnOldAllele = max_MOI * i + oldalleles.shape[1]
            endColumnNewAllele = max_MOI * i + newalleles.shape[1]
            self.alleles0[:, startColumn:endColumnOldAllele] = oldalleles[
                genotypedata_RR["Sample ID"].str.contains("Day 0"), :
            ]
            self.allelesf[:, startColumn:endColumnOldAllele] = oldalleles[
                genotypedata_RR["Sample ID"].str.contains("Day Failure"), :
            ]
            self.recoded0[:, startColumn:endColumnNewAllele] = newalleles[
                genotypedata_RR["Sample ID"].str.contains("Day 0"), :
            ]
            self.recodedf[:, startColumn:endColumnNewAllele] = newalleles[
                genotypedata_RR["Sample ID"].str.contains("Day Failure"), :
            ]

    @classmethod
    def recode_additional_neutral(
        cls,
        additional_neutral: pd.DataFrame,
        alleles_definitions_RR: pd.DataFrame,
        locinames: np.ndarray,
        max_MOI: int):
        '''
        TODO: Verify what this actually does?
        If additional_neutral data is present, return a recoded version
        TODO: Figure out what recoding actually means?

        :param additional_neutral: The genotype data for this site the sample
        data
        :param alleles_definitions_RR: TODO:
        :param locinames: The names of the loci in additional_neutral to get
        alleles for
        :param max_MOI: The maximum multiplicity of infection in the additional
        data
        '''
        if additional_neutral.size == 0 or additional_neutral.shape[0] == 0:
            return None

        recoded_additional_neutral = np.zeros((additional_neutral.shape[0], max_MOI * locinames.size))
        for i, locus in enumerate(locinames):
            oldalleles, newalleles = cls._get_original_alleles(
                additional_neutral, alleles_definitions_RR, locus, i)

            startColumn = max_MOI * (i - 1)  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
            endColumnOldAllele = max_MOI * (i - 1) + oldalleles.shape[1]
            recoded_additional_neutral[:, startColumn:endColumnOldAllele] = newalleles
        return recoded_additional_neutral

    @classmethod
    def _get_original_alleles(
        cls,
        samples_df: pd.DataFrame,
        alleles_definitions_RR: pd.DataFrame,
        locus: str,
        locus_index: int):
        '''
        TODO: Verify what this actually does?
        Return the oldalleles/newalleles from the given dataframe, with NaN
        values set to 0

        :param samples_df: The dataframe containing the samples
        :param alleles_definitions_RR:
        :param locus: The locus in the dataframe to get alleles for
        :param locus_index: The index of the locus among the unique locinames in
        the dataframe
        :return: A tuple of 2 numpy arrays: (oldalleles, newalleles)
        '''
        locicolumns = samples_df.columns.str.contains(f"{locus}_")

        oldalleles = samples_df.loc[:, locicolumns].to_numpy()
        newalleles = np.copy(oldalleles)
        num_columns = oldalleles.shape[1]
        for j in range(num_columns):
            newalleles[:, j] = np.array(list(map(
                lambda x: cls._recode_allele(
                    alleles_definitions_RR[locus_index].to_numpy(),
                    oldalleles[x, j]),
                range(0, oldalleles.shape[0])
            )))

        # Set all nans in either array to 0
        oldalleles[np.isnan(oldalleles)] = 0
        oldalleles[np.isnan(newalleles)] = 0
        newalleles[np.isnan(newalleles)] = 0

        return oldalleles, newalleles

    @classmethod
    def _recode_allele(
        cls,
        alleles_definitions_subset: np.ndarray,
        proposed: np.ndarray):
        '''
        TODO: Figure out what this actually does?

        :param alleles_definitions_subset: 2D numpy array (?)
        :param proposed: 1D numpy vector (?)
        :return: Returns a single integer index, or np.nan if no valid index was
        found
        '''
        # verify shapes.
        if len(alleles_definitions_subset.shape) != 2 or alleles_definitions_subset.shape[1] != 2:
            raise ValueError(f'Improper alleles_definition_subset shape {alleles_definitions_subset.shape} (expected (:,2))')
        if len(proposed.shape) > 1:
            raise ValueError(f'Improper proposed vector shape {proposed.shape} (expected 1D vector)')

        # alleles_definitions_subset ranges guaranteed to be non-overlapping, so it will always fall within either 0 or exactly 1 of the ranges (i.e. rows)
        result = np.argwhere(np.logical_and(proposed > alleles_definitions_subset[:, 0], proposed <= alleles_definitions_subset[:, 1]))
        result = result.reshape(-1)
        if (result.size == 0):
            result = np.nan
        else:
            return result[0]
        return result

    def randomize_initial_assignments(
        self,
        num_ids: int,
        num_loci: int,
        max_MOI: int,
        alleles_definitions_RR: pd.DataFrame,
        seed: int = None):
        '''
        TODO: Elaborate
        Assign random initial values for the hidden alleles and
        reinfection/recrudescence classifications, based on the prior calculated
        frequencies for the malaria dataset
        '''
        random = np.random.RandomState(seed)
        for id_index in range(num_ids):
            # 50% chance if this sample should be initialized as a reinfection
            # or recrudescence
            if random.uniform(size=1):
                self.classification[id_index] = SampleClass.RECRUDESCENCE
            for locus_index in range(num_loci):
                self._randomize_hidden_alleles(
                    id_index,
                    locus_index,
                    max_MOI,
                    alleles_definitions_RR,
                    random)
                self._randomize_recrudescences(id_index, locus_index, max_MOI)

    def _randomize_hidden_alleles(
        self,
        id_index: int,
        locus_index: int,
        max_MOI: int,
        alleles_definitions_RR: pd.DataFrame,
        rand: np.random.RandomState):
        '''
        TODO: Elaborate
        Randomizes the initial allele values for the given ID/locus

        :param id_index: The index of the current ID to set alleles for
        :param locus_index: The index of the current locus to set alleles for
        :param max_MOI: The maximum multiplicity of infection in the dataset
        :param alleles_definitions_RR:
        :param rand: A generator for the random numbers in this function
        '''
        self._randomize_allele_group(
            self.alleles0, self.recoded0, self.hidden0, self.MOI0,
            id_index, locus_index, max_MOI, alleles_definitions_RR, rand)
        self._randomize_allele_group(
            self.allelesf, self.recodedf, self.hiddenf, self.MOIf,
            id_index, locus_index, max_MOI, alleles_definitions_RR, rand)

    def _randomize_allele_group(
        self,
        alleles: np.ndarray,
        recoded: np.ndarray,
        hidden: np.ndarray,
        MOIs: np.ndarray,
        id_index: int,
        locus_index: int,
        max_MOI: int,
        alleles_definitions_RR: pd.DataFrame,
        rand: np.random.RandomState):
        '''
        TODO: Cut down on the parameter list (combine last few elements in
        tuple?)
        NOTE: Depends on assuming alleles, recoded, and hidden will be modified
        in-place (i.e. that they're passed by reference)
        '''
        i = id_index
        j = locus_index
        # TODO: Start/end of what? The portion of the row w/ this locus information?
        start = max_MOI * j
        end = max_MOI * (j + 1)

        num_alleles = np.count_nonzero(alleles[i, start:end])
        num_missing = MOIs[i] - num_alleles

        missing_alleles_indices = np.arange(start, end)[
            np.where(alleles[i, start: start + MOIs[i]] == 0)
        ]
        present_alleles_indices = np.delete(np.arange(start, end), missing_alleles_indices)

        # Sample to randomly initialize the alleles/hidden variables
        if num_alleles > 0:
            hidden[i, present_alleles_indices] = 0
        elif num_missing > 0:
            new_hidden_alleles = rand.choice(
                np.arange(
                    0, int(self.frequencies_RR[0][j])
                ),  # Select from first row (count of how many probabilities they are)
                size=num_missing,
                replace=True,
                p=self.frequencies_RR[1][j, 0: int(self.frequencies_RR[0][j])]
            )  # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
            recoded[i, missing_alleles_indices] = new_hidden_alleles
            # calculate row means
            alleles[i, missing_alleles_indices] = np.mean(alleles_definitions_RR[j], axis=1)[
                new_hidden_alleles
            ]  # hidden alleles get mean allele length
            hidden[i, missing_alleles_indices] = 1

    def _randomize_recrudescences(
        self,
        id_index: int,
        locus_index: int,
        max_MOI: int):
        '''
        TODO: Elaborate (since this doesn't seem to actually do anything random?
        So what is it doing?)
        Randomly assign initial recrudscences/reinfections
        '''
        i = id_index
        j = locus_index

        # determine which alleles are recrudescing (for beginning, choose closest pair)
        allpossiblerecrud = np.stack(
            np.meshgrid(np.arange(self.MOI0[i]), np.arange(self.MOIf[i]))
        ).T.reshape(-1, 2)

        # TODO: Much of the below code near-duplicated?
        allele0_col_indices = max_MOI * j + allpossiblerecrud[:, 0]
        allelef_col_indices = max_MOI * j + allpossiblerecrud[:, 1]

        recrud_distances = np.abs(
            self.alleles0[i, allele0_col_indices] - self.allelesf[i, allelef_col_indices]
        )

        closest_recrud_index = np.argmin(recrud_distances)

        self.mindistance[i, j] = recrud_distances[closest_recrud_index]
        self.alldistance[i, j, :recrud_distances.size] = recrud_distances

        self.allrecrf[i, j, :allpossiblerecrud.shape[0]] = self.recodedf[
            i, max_MOI * j + allpossiblerecrud[:, 1]
        ]

        def set_recrudescences(recr, recr_repeats, recoded, is_0=True):
            # TODO: Verify what the purpose of this actually is?
            recrud_column = 0 if is_0 else 1
            recr[i, j] = max_MOI * j + allpossiblerecrud[closest_recrud_index, recrud_column]
            recr_repeats[i, j] = np.sum(
                recoded[i, max_MOI * j: max_MOI * (j+1)] == recoded[i, int(recr[i, j])])

        set_recrudescences(self.recr0, self.recr_repeats0, self.recoded0)
        set_recrudescences(self.recrf, self.recr_repeatsf, self.recodedf, False)
