from collections import namedtuple
from typing import List

import bottleneck as bn
import numpy as np
import pandas as pd

from api.site_instance_state import HiddenAlleleType, SampleType, SiteInstanceState


# State to hold shared local variables across calculations
SwitchHiddenState = namedtuple('_SwitchState', 'z, chosen, chosenlocus, is_chosen_valid, old, new, newallele_length, oldalleles, allpossiblerecrud')


def switch_hidden(x, nloci, maxMOI, alleles_definitions_RR_arr: List[np.ndarray], state: SiteInstanceState, rand: np.random.RandomState):
    '''
    Update the hidden/inferred alleles for the samples

    NOTE: Update takes place entirely via side-effects/in-place updates to the
    state (this function doesn't currently return anything)

    :param x: The index of the sample ID to update
    :param nloci: The number of loci in the dataset
    :param maxMOI: The maximum multiplicity of infection in the whole dataset
    :param alleles_definitions_RR_arr: ?
    :param state: The current state of the algorithm
    :param rand: The random number generator to use
    '''

    # Section A: If number of inferred alleles > 0
    # It will probably be more efficient to sum the two seperately, because concatenation
    # could induce memory-related performance cost, if a new memory block is being created behind the scenes.
    inferred_allele_count = bn.nansum(state.hidden0) + bn.nansum(state.hiddenf)
    if inferred_allele_count <= 0:
        return

    sh_state = _setup_initial_state(x, nloci, maxMOI, alleles_definitions_RR_arr, state, rand)
    is_reinfection = state.classification[x] == SampleType.REINFECTION.value
    if is_reinfection:
        _update_reinfection(state, sh_state, x, maxMOI)
    else:
        _update_recrudescence(state, sh_state, x, maxMOI)


def _setup_initial_state(x, nloci, maxMOI, alleles_definitions_RR_arr: List[np.ndarray], state: SiteInstanceState, rand: np.random.RandomState) -> SwitchHiddenState:
    '''
    Sets up the initial local variables switch_hidden uses throughout its calculations, including the following:
        z - A random probability between 0 and 1 that's used to determine if we should update the allele
        chosen - The randomly-selected index of the allele to update
        chosenlocus - The locus of the allele we're updating
        is_chosen_valid - A boolean indicating if the chosen allele index fits within the state's allele arrays
        old - ?
        new - ?
        newallele_length - The length to update the allele to (?)
        oldalleles - The bin indices of all the hidden alleles
        allpossiblerecrud - An array of all possible combinations of recrudescing alleles

    :param x: The index of the sample ID to update
    :param nloci: The number of loci in the dataset
    :param maxMOI: The maximum multiplicity of infection in the whole dataset
    :param alleles_definitions_RR_arr: ?
    :param state: The current state of the algorithm
    :param rand: The random number generator to use
    :return: An initialized SwitchHiddenState tuple
    '''
    z = rand.uniform(size=1)

    inferred_allele_indices = np.where(np.concatenate((state.hidden0[x], state.hiddenf[x])) == HiddenAlleleType.MISSING.value)[0]

    # Figure out chosen & chosenlocus
    # QUESTION: What is 'chosen'? Chosen..... allele?
    chosen = rand.choice(inferred_allele_indices) + 1

    # TODO: Rename this. what is is_chosen_valid?
    is_chosen_valid = chosen <= (nloci * maxMOI)

    if not is_chosen_valid:
        chosen -= nloci * maxMOI
    chosenlocus = int(np.ceil(chosen / maxMOI))

    # Subtract by 1 to index from 0
    chosen -= 1
    chosenlocus -= 1

    old = state.recoded0[x, chosen].astype(np.int64)
    new = rand.choice(np.arange(state.frequencies_RR.lengths[chosenlocus])) + 1
    if is_chosen_valid:
        oldalleles = _get_old_alleles(state.recoded0, state.hidden0, x, chosenlocus, maxMOI)
    else:
        oldalleles = _get_old_alleles(state.recodedf, state.hiddenf, x, chosenlocus, maxMOI)

    newallele_length = np.mean(alleles_definitions_RR_arr[chosenlocus][new-1, :]) + rand.normal(loc=0, scale=state.frequencies_RR.variability[chosenlocus], size=1)

    allpossiblerecrud = state.get_all_possible_recrud(sample_id=x)

    sh_state = SwitchHiddenState(z, chosen, chosenlocus, is_chosen_valid, old, new, newallele_length, oldalleles, allpossiblerecrud)
    return sh_state


# TODO: Refactor out commonalities between reinfection/recrudescence update
def _update_reinfection(state: SiteInstanceState, sh: SwitchHiddenState, x: int, maxMOI: int):
    '''
    Handle updating the hidden variables for a reinfected sample

    :param state: The current state of the algorithm
    :param sh: The computed local variables within switch hidden to operate on
    :param x: The sample ID to update
    :param maxMOI: The maximum multiplicity of infection for the overall dataset
    '''
    repeatednew = 1 if np.any(sh.oldalleles == sh.new) else state.qq

    numerator = np.sum(state.frequencies_RR.matrix[sh.chosenlocus, 0:state.frequencies_RR.lengths[sh.chosenlocus]] * state.dvect[state.correction_distance_matrix[sh.chosenlocus][sh.new - 1].astype(np.int64)]) * repeatednew
    denominator = np.sum(state.frequencies_RR.matrix[sh.chosenlocus, 0:state.frequencies_RR.lengths[sh.chosenlocus]] * state.dvect[state.correction_distance_matrix[sh.chosenlocus][sh.old].astype(np.int64)]) * repeatednew

    alpha = numerator / denominator if denominator != 0 else 0
    if sh.z >= alpha:
        return

    if (sh.is_chosen_valid):
        state.recoded0[x, sh.chosen] = sh.new - 1
    else:
        state.recodedf[x, sh.chosen] = sh.new - 1
    state.alleles0[x, sh.chosen] = sh.newallele_length

    # Not sure what this list actually is, I suspect it refers to distances between
    # microsatelites of alleles. Either way, we calculate it twice, and the calculation seems
    # kinda expensive, so I'm defining it here to avoid that issue.
    dist_list = list(map(lambda y: _unknownhelper_1(state, x, maxMOI, sh.chosenlocus, sh.allpossiblerecrud, y), np.arange(0, sh.allpossiblerecrud.shape[0])))

    closestrecrud = np.argmin(dist_list)
    state.mindistance[x, sh.chosenlocus] = dist_list[closestrecrud]
    state.alldistance[x, sh.chosenlocus, :sh.allpossiblerecrud.shape[0]] = dist_list
    state.allrecrf[x, sh.chosenlocus, :sh.allpossiblerecrud.shape[0]] = state.recodedf[x, maxMOI * sh.chosenlocus + sh.allpossiblerecrud[:, 1]]
    state.recr0[x, sh.chosenlocus] = maxMOI * sh.chosenlocus + sh.allpossiblerecrud[:, 0][closestrecrud]
    state.recrf[x, sh.chosenlocus] = maxMOI * sh.chosenlocus + sh.allpossiblerecrud[:, 1][closestrecrud]


def _update_recrudescence(state: SiteInstanceState, sh: SwitchHiddenState, x: int, maxMOI: int):
    '''
    Handle updating the hidden variables for a recrudescing sample

    :param state: The current state of the algorithm
    :param sh: The computed local variables within switch hidden to operate on
    :param x: The sample ID to update
    :param maxMOI: The maximum multiplicity of infection for the overall dataset
    '''
    repeatedold = 1 if np.any(sh.oldalleles == sh.old) else state.qq
    repeatednew = 1 if np.any(sh.oldalleles == sh.new) else state.qq

    newclosestrecrud, newmindistance, newalldistance, newallrecrf = _calculate_new_distances(state, sh, x, maxMOI)

    ## likelihoodnew
    likelihoodnew_numerator = state.dvect[np.round(newalldistance).astype(np.int64)]
    likelihoodnew_demominator = np.sum(state.frequencies_RR.matrix[sh.chosenlocus][:state.frequencies_RR.lengths[sh.chosenlocus]] * state.dvect[state.correction_distance_matrix[sh.chosenlocus][newallrecrf[0].astype(np.int64)].astype(np.int64)])
    likelihoodnew = np.nanmean(likelihoodnew_numerator / likelihoodnew_demominator) * repeatednew

    ## likelihoodold
    temp = np.round(state.alldistance[x, sh.chosenlocus])
    temp = temp[~np.isnan(temp)].astype(np.int64)
    likelihoodold_numerator = state.dvect[temp]

    ### NEED TO REVISIT
    temp_allrecrf = state.allrecrf[x, sh.chosenlocus, :maxMOI**2]
    temp_allrecrf = temp_allrecrf[~np.isnan(temp_allrecrf)].astype(np.int64)

    temp = state.frequencies_RR.matrix[sh.chosenlocus, :state.frequencies_RR.lengths[sh.chosenlocus]] * state.dvect[state.correction_distance_matrix[sh.chosenlocus][temp_allrecrf].astype(np.int64)]
    likelihoodold_denominator = []
    for i in temp:
        likelihoodold_denominator.append(np.sum(i))

    likelihoodold_denominator = np.asarray(likelihoodold_denominator)
    likelihoodold = np.nanmean(likelihoodold_numerator / likelihoodold_denominator) * repeatedold

    # Prepare an 'alpha' from our previous likelihoods that will determine whether or not
    # We are going to update the state.
    alpha = 0
    if likelihoodnew == likelihoodold:
        alpha = 1
    elif not likelihoodold == 0:
        alpha = likelihoodnew / likelihoodold

    if sh.z >= alpha:
        return

    state.recoded0[x, sh.chosen] = sh.new - 1
    if sh.is_chosen_valid:
        state.alleles0[x, sh.chosen] = sh.newallele_length
    else:
        state.allelesf[x, sh.chosen] = sh.newallele_length
    state.mindistance[x, sh.chosenlocus] = newmindistance
    state.alldistance[x, sh.chosenlocus, :sh.allpossiblerecrud.shape[0]] = newalldistance
    state.allrecrf[x, sh.chosenlocus, :sh.allpossiblerecrud.shape[0]] = newallrecrf
    state.recr0[x, sh.chosenlocus] = maxMOI * (sh.chosenlocus) + sh.allpossiblerecrud[:, 0][newclosestrecrud]
    state.recrf[x, sh.chosenlocus] = maxMOI * (sh.chosenlocus) + sh.allpossiblerecrud[:, 1][newclosestrecrud]


def _get_old_alleles(recoded: np.ndarray, hidden: np.ndarray, id_index: int, chosen_locus: int, max_MOI: int) -> np.ndarray:
    '''
    Returns the bin indices (in recoded) of the previous inferred alleles

    :param recoded: The range bin indices of each allele
    :param hidden: An array classifying each recoded allele as "hidden" or "observed"
    :param id_index: The sample ID number to get alleles for
    :param chosen_locus: The locus index to get alleles for
    :param max_MOI: The max multiplicity of infection in the dataset
    '''
    start_index = chosen_locus * max_MOI
    end_index = start_index + max_MOI
    return recoded[id_index, start_index: end_index][
        hidden[id_index, start_index: end_index] == HiddenAlleleType.MISSING.value
    ]

def _calculate_new_distances(state: SiteInstanceState, sh: SwitchHiddenState, id_index: int, max_MOI: int):
    '''
    Returns the updated distance-based parameters for state
    '''
    if sh.is_chosen_valid:
        tempalleles = state.alleles0[id_index, max_MOI * sh.chosenlocus: max_MOI * (sh.chosenlocus + 1)]
        tempalleles[sh.chosen - sh.chosenlocus * max_MOI] = sh.newallele_length
        temprecoded = state.recoded0[id_index, max_MOI * sh.chosenlocus: max_MOI * (sh.chosenlocus + 1)]
        temprecoded[sh.chosen - sh.chosenlocus * max_MOI] = sh.new - 1

        # This was seen in the earlier branch as well, I believe
        # this is a distance calculation taking into account that we have
        # no day0 data to work with. Either way, it was calculated twice,
        # so I have put the list here again.

        newalldistance = list(map(lambda y: _unknownhelper_2(state,tempalleles,id_index,max_MOI,sh.chosenlocus,sh.allpossiblerecrud,y), np.arange(0, sh.allpossiblerecrud.shape[0])))

        newclosestrecrud = np.argmin(newalldistance)
        newmindistance = newalldistance[newclosestrecrud]
        newallrecrf = state.recodedf[id_index, max_MOI * sh.chosenlocus + sh.allpossiblerecrud[1]]
    else:
        tempalleles = state.allelesf[id_index, max_MOI * sh.chosenlocus: max_MOI * (sh.chosenlocus + 1)]
        tempalleles[sh.chosen - sh.chosenlocus * max_MOI] = sh.newallele_length
        temprecoded = state.recodedf[id_index, max_MOI * sh.chosenlocus: max_MOI * (sh.chosenlocus + 1)]
        temprecoded[sh.chosen - sh.chosen_locus * max_MOI] = sh.new - 1

        newalldistance = list(map(lambda y:_unknownhelper_3(state,tempalleles,id_index,max_MOI,sh.chosenlocus,sh.allpossiblerecrud,y) , np.arange(0, sh.allpossiblerecrud.shape[0])))

        newclosestrecrud = np.argmin(newalldistance)
        newmindistance = newalldistance[newclosestrecrud]
        newallrecrf = temprecoded[sh.allpossiblerecrud[:, 1]]

    return newclosestrecrud, newmindistance, newalldistance, newallrecrf


# TODO: Determine what this subroutine actually is.
# I think its for calculating a distance between microsatelite alleles. Is that right?
def _unknownhelper_1(state,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[y, 0]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[y, 1]])
    return value


# TODO: Determine what this subroutine actually is.
# I think its for calculating a distance between microsatelite alleles, when we don't have day 0 data. Is that right?
def _unknownhelper_2(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(tempalleles[allpossiblerecrud[y, 0]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[y, 1]])
    return value


# TODO: Determine what this subroutine actually is.
# Looks similar to our other evaluations, but we use day0 data. I guess this is the case where we don't have f data, but we have 0.
def _unknownhelper_3(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(tempalleles[allpossiblerecrud[y, 1]] - state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[y, 0]])
    return value
