import enum
from typing import List

import numpy as np
import pandas as pd

from api.calculate_frequencies import calculate_frequencies3


class SampleType(enum.Enum):
    REINFECTION = 0
    RECRUDESCENCE = 1


class HiddenAlleleType(enum.Enum):
    # TODO: Verify I don't have these switched?
    MISSING = 0
    OBSERVED = 1
    UNKNOWN = np.nan    # TODO: Confirm this is correct?


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
        alleles_definitions_RR: List[pd.DataFrame]):
        '''
        Sets up the deterministic initial state of the algorithm
        TODO: Elaborate, shorten argument list?
        '''
        self._create_empty_state(ids.size, locinames.size, maxMOI)

        self.MOI0, self.MOIf = self._calculate_sample_MOI(
            genotypedata_RR, ids, locinames)
        self._initialize_alleles(
            genotypedata_RR, alleles_definitions_RR, locinames, maxMOI)
        self.recoded_additional_neutral = self.recode_additional_neutral(
            additional_neutral, alleles_definitions_RR, locinames, maxMOI)
        self.dvect = self._get_initial_dvect(alleles_definitions_RR)
        self.qq = np.nan
        self.dposterior = 0.75  # TODO: What is this?

        # TODO: Unsure if these should be here, since they're read-only state?
        # estimate frequencies
        self.frequencies_RR = calculate_frequencies3(
            pd.concat([genotypedata_RR, additional_neutral]),
            alleles_definitions_RR)
        self.correction_distance_matrix = self._get_correction_distances(
            alleles_definitions_RR)

    def _create_empty_state(self, num_ids: int, num_loci: int, max_MOI: int):
        '''
        Creates the initial empty data structures in the state to be populated
        later. The data structures are as follows:

        alleles0 - Holds the allele fragment lengths found on Day 0
        recoded0 - Holds the indices of the bins the allele0 lengths fall into
        hidden0 - Records if a given allele was directly observed in the
        original dataset (0) or has to be inferred (1) (TODO: What does nan
        mean?)
        recr0 - Records which strain we think is recrudescing for each sample
        (e.g. recr0[i, j] = 2 means "For sample i at locus j, we think the
        2nd strain of the locus is the recrudescent one")
        *f - All the "f" variants are the same thing, but for the day of failure
        mindistance - The closest of all possible allele length pairs for a
        given sample/locus
        alldistance - All possible allele length pairs for a given sample/locus
        allrecrf - Holds all possible allele combinations for the given sample/
        locus, and marks which one we estimate is the recrudescent one
        classification - Holds whether we think a given sample is a reinfection
        or a recrudescence

        :param num_ids: The number of samples in the dataset
        :param num_loci: The number of loci in the dataset
        :param max_MOI: The max multiplicity of infection in the dataset
        '''
        self.alleles0 = np.zeros((num_ids, max_MOI * num_loci))
        self.recoded0 = np.zeros((num_ids, max_MOI * num_loci))
        self.hidden0 = np.full_like(np.empty((num_ids, max_MOI * num_loci)),
            HiddenAlleleType.UNKNOWN.value)
        self.recr0 = np.full_like(np.empty((num_ids, num_loci)), np.nan)
        self.allelesf = np.copy(self.alleles0)
        self.recodedf = np.copy(self.recoded0)
        self.hiddenf = np.copy(self.hidden0)
        self.recrf = np.copy(self.recr0)
        self.mindistance = np.zeros((num_ids, num_loci))
        self.alldistance = np.full_like(np.empty((num_ids, num_loci, max_MOI ** 2)), np.nan)
        self.allrecrf = np.full_like(np.empty((num_ids, num_loci, max_MOI ** 2)), np.nan)
        self.classification = np.repeat(SampleType.REINFECTION.value, num_ids)

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

        :param genotypedata_RR: The day 0 and day of failure genotype data for
        this site's samples
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
                        genotypedata_RR["Sample ID"].str.contains(f"{ID} Day 0"), locicolumns
                    ].isna()
                )
                num_allelesf = np.count_nonzero(
                    ~genotypedata_RR.loc[
                        genotypedata_RR["Sample ID"].str.contains(f"{ID} Day Failure"),
                        locicolumns,
                    ].isna()
                )

                MOI0[i] = np.max([MOI0[i], num_alleles0])
                MOIf[i] = np.max([MOIf[i], num_allelesf])
        return MOI0, MOIf

    def _initialize_alleles(
        self,
        genotypedata_RR: pd.DataFrame,
        alleles_definitions_RR: List[pd.DataFrame],
        locinames: np.ndarray,
        max_MOI: int):
        '''
        TODO: Verify what this actually does?
        Initialize the alleles/recoded state, filling them with the appropriate
        initial/failure data from the dataframe

        :param genotypedata_RR: The day 0 and day of failure genotype data for
        this site's samples
        :param alleles_definitions_RR: The bins/ranges that the alleles lengths
        could fall into
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
        alleles_definitions_RR: List[pd.DataFrame],
        locinames: np.ndarray,
        max_MOI: int):
        '''
        TODO: Verify what this actually does?
        If additional_neutral data is present, return a recoded version with
        the allele lengths appropriately binned; otherwise, returns an empty
        numpy array

        :param additional_neutral: The background genotype data for this site's
        sample data
        :param alleles_definitions_RR: TODO:
        :param locinames: The names of the loci in additional_neutral to get
        alleles for
        :param max_MOI: The maximum multiplicity of infection in the additional
        data
        :return: The recoded additional_neutral array, or an empty numpy array
        if there is no additional_neutral data
        '''
        if additional_neutral.size == 0 or additional_neutral.shape[0] == 0:
            return np.empty([0,0])

        recoded_additional_neutral = np.zeros((additional_neutral.shape[0], max_MOI * locinames.size))
        for i, locus in enumerate(locinames):
            oldalleles, newalleles = cls._get_original_alleles(
                additional_neutral, alleles_definitions_RR, locus, i)

            # TODO: Same indexing as for _initializing_alleles?
            startColumn = max_MOI * i
            endColumn = startColumn + oldalleles.shape[1]
            recoded_additional_neutral[:, startColumn:endColumn] = newalleles
        return recoded_additional_neutral

    @classmethod
    def _get_original_alleles(
        cls,
        samples_df: pd.DataFrame,
        alleles_definitions_RR: List[pd.DataFrame],
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
        proposed: float):
        '''
        Returns the index of the alleles_definitions bin that the proposed
        allele length falls into, or np.nan if it doesn't fall within any
        of the subset's ranges

        :param alleles_definitions_subset: Nx2 2D numpy array, representing
        (mutually exclusive) ranges of allele lengths the proposed value can
        fall between
        :param proposed: Value that should fall between one of the ranges in the
        subset
        :return: Returns a single integer index of the row/range the proposed
        number falls within in the subset, or np.nan if no valid range was
        found
        '''
        # verify shapes.
        if len(alleles_definitions_subset.shape) != 2 or alleles_definitions_subset.shape[1] != 2:
            raise ValueError(f'Improper alleles_definition_subset shape {alleles_definitions_subset.shape} (expected (:,2))')

        # alleles_definitions_subset ranges guaranteed to be non-overlapping, so it will always fall within either 0 or exactly 1 of the ranges (i.e. rows)
        result = np.argwhere(np.logical_and(proposed > alleles_definitions_subset[:, 0], proposed <= alleles_definitions_subset[:, 1]))
        result = result.reshape(-1)
        if (result.size == 0):
            result = np.nan
        else:
            return result[0]
        return result

    @classmethod
    def _get_initial_dvect(cls, alleles_definitions_RR: pd.DataFrame):
        '''
        TODO: Understand this better?
        Return the initial distance vector (estimating the likelihood of error
        in the analysis)
        '''
        ranges = []
        for dataframe in alleles_definitions_RR:
            # Get the range (max-min) of the first "nloci" dataframes, then the max of all those
            ranges.append(dataframe.max().max() - dataframe.min().min())

        dvect = np.zeros(1 + int(round(max(ranges))))
        dvect[0] = 0.75
        dvect[1] = 0.2
        dvect[2] = 0.05

        return dvect

    @classmethod
    def _get_correction_distances(
        cls,
        alleles_definitions_RR: List[pd.DataFrame]) -> List[np.ndarray]:
        '''
        TODO: Verify this description
        Returns the matrix of distances between each pair of alleles, for each
        locus
        :param alleles_definitions_RR: TODO:
        :return: Python array of 2D matrices of distances between each allele
        pair for a given locus
        '''
        correction_distance_matrix = [] # for each locus, matrix of distances between each allele
        # TODO: Vectorize this (it seems fairly doable)
        for i in range(len(alleles_definitions_RR)):
            # Wrap mean call in "array" so we get a 2D array we can transpose (getting us a grid of distances, not just a 1D vector)
            distances = np.array([np.mean(alleles_definitions_RR[i], axis=1)])
            distance_combinations = np.abs(distances.T - distances)
            correction_distance_matrix.append(distance_combinations)
        return correction_distance_matrix

    def randomize_initial_assignments(
        self,
        num_ids: int,
        num_loci: int,
        max_MOI: int,
        alleles_definitions_RR: List[pd.DataFrame],
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
                self.classification[id_index] = SampleType.RECRUDESCENCE.value
            for locus_index in range(num_loci):
                self._randomize_hidden_alleles(
                    id_index,
                    locus_index,
                    max_MOI,
                    alleles_definitions_RR,
                    random)
                self._assign_closest_recrudescences(id_index, locus_index, max_MOI)

        self.qq = self._get_initial_qq(self.hidden0, self.hiddenf)

    def _randomize_hidden_alleles(
        self,
        id_index: int,
        locus_index: int,
        max_MOI: int,
        alleles_definitions_RR: List[pd.DataFrame],
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
        alleles_definitions_RR: List[pd.DataFrame],
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
            np.where(alleles[i, start: start + MOIs[i]] == HiddenAlleleType.MISSING.value)
        ]
        present_alleles_indices = np.delete(np.arange(start, end), missing_alleles_indices)

        # Sample to randomly initialize the alleles/hidden variables
        if num_alleles > 0:
            hidden[i, present_alleles_indices] = HiddenAlleleType.MISSING.value
        if num_missing == 0:
            return

        new_hidden_alleles = rand.choice(
            # Select from first row (count of how many probabilities they are)
            np.arange(0, int(self.frequencies_RR[0][j])),
            size=num_missing,
            replace=True,
            p=self.frequencies_RR[1][j, 0: int(self.frequencies_RR[0][j])]
        )
        # Choose random initial data for missing alleles
        recoded[i, missing_alleles_indices] = new_hidden_alleles
        # calculate row means (mean allele lengths)
        alleles[i, missing_alleles_indices] = np.mean(alleles_definitions_RR[j], axis=1)[new_hidden_alleles]
        hidden[i, missing_alleles_indices] = HiddenAlleleType.OBSERVED.value

    def _assign_closest_recrudescences(
        self,
        id_index: int,
        locus_index: int,
        max_MOI: int):
        '''
        Finds the closest possible recrudescing allele pairs known for each
        locus and sample ID, records them as our initial guess for the
        recrudescence, and updates the variables appropriately.

        The closest allele pair is assigned as our initial guess for the most
        likely allele to be recrudescent with the sample's day 0 allele.

        :param id_index: The index of the sample ID being evaluated
        :param locus_index: The index of the locus being evaluated
        :max_MOI: The maximum multiplicity of infection for the dataset
        '''
        i = id_index
        j = locus_index

        # calculate all possible pairs of alleles that could be recrudescing for
        # this sample
        # NOTE: Correct indices generated, but in a different order than R code
        allpossiblerecrud = np.stack(
            np.meshgrid(np.arange(self.MOI0[i]), np.arange(self.MOIf[i]))
        ).T.reshape(-1, 2)

        # TODO: Much of the below code near-duplicated?
        allele0_col_indices = max_MOI * j + allpossiblerecrud[:, 0]
        allelef_col_indices = max_MOI * j + allpossiblerecrud[:, 1]

        # calculate distances between each possible pair of alleles
        recrud_distances = np.abs(
            self.alleles0[i, allele0_col_indices] - self.allelesf[i, allelef_col_indices]
        )

        # select the closest pair, and record the closest distance/distances for
        # this sample
        closest_recrud_index = np.argmin(recrud_distances)

        self.mindistance[i, j] = recrud_distances[closest_recrud_index]
        self.alldistance[i, j, :recrud_distances.size] = recrud_distances

        self.allrecrf[i, j, :allpossiblerecrud.shape[0]] = self.recodedf[
            i, max_MOI * j + allpossiblerecrud[:, 1]
        ]

        def set_recrudescences(recr, is_day_0=True):
            # TODO: Verify what the purpose of this actually is?
            recrud_column = 0 if is_day_0 else 1
            recr[i, j] = max_MOI * j + allpossiblerecrud[closest_recrud_index, recrud_column]

        set_recrudescences(self.recr0)
        set_recrudescences(self.recrf, False)

    @classmethod
    def _get_initial_qq(cls, hidden0: np.ndarray, hiddenf: np.ndarray):
        '''
        TODO: What does qq stand for?
        Initial estimate of q, the probability of an allele being missed

        :param hidden0: TODO:
        :param hiddenf: TODO:
        :return: A single number q, the mean of the known hidden variables
        '''
        return np.nanmean(np.concatenate([hidden0, hiddenf]))
