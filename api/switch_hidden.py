import itertools

import numpy as np
import pandas as pd

from api.site_instance_state import HiddenAlleleType, SampleType

def switch_hidden(x, nloci, maxMOI, alleles_definitions_RR, state, rand: np.random.RandomState):
    z = rand.uniform(size=1)

    # Section A: If number of inferred alleles > 0
    # It will probably be more efficient to sum the two seperately, because concatenation
    # could induce memory-related performance cost, if a new memory block is being created behind the scenes.
    inferred_allele_count = np.nansum(state.hidden0) + np.nansum(state.hiddenf)
    if inferred_allele_count <= 0:
        return

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

    is_reinfection = state.classification[x] == SampleType.REINFECTION.value
    if is_reinfection:
        old = state.recoded0[x, chosen].astype(np.int64)
        new = rand.choice(np.arange(state.frequencies_RR[0][chosenlocus])) + 1
        oldalleles = state.recoded0[x, np.intersect1d(np.arange((chosenlocus * maxMOI), (chosenlocus + 1) * maxMOI), np.where(state.hidden0[x] == HiddenAlleleType.MISSING.value)[0])]
        repeatednew = state.qq

        if np.sum(oldalleles == new) >= 1:
            repeatednew = 1

        numerator = np.sum(state.frequencies_RR[1][chosenlocus, 0:state.frequencies_RR[0][chosenlocus]] * state.dvect[state.correction_distance_matrix[chosenlocus][new - 1].astype(np.int64)]) * repeatednew
        denominator = np.sum(state.frequencies_RR[1][chosenlocus, 0:state.frequencies_RR[0][chosenlocus]] * state.dvect[state.correction_distance_matrix[chosenlocus][old].astype(np.int64)]) * repeatednew

        alpha = numerator / denominator if denominator != 0 else 0
        if z >= alpha:
            return

        if (is_chosen_valid):
            state.recoded0[x, chosen] = new - 1
        else:
            state.recodedf[x, chosen] = new - 1
        newallele_length = (np.mean((alleles_definitions_RR[chosenlocus]["0"][new-1], alleles_definitions_RR[chosenlocus]["1"][new-1])) + rand.normal(0, state.frequencies_RR[2][chosenlocus], 1))[0]
        state.alleles0[x, chosen] = newallele_length

        inputVectors = list(itertools.product(np.arange(state.MOIf[x]), np.arange(state.MOI0[x])))
        allpossiblerecrud = pd.DataFrame(inputVectors)
        order = [1, 0] # setting column's order
        allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
        allpossiblerecrud.columns = [0, 1]

        # Not sure what this list actually is, I suspect it refers to distances between
        # microsatelites of alleles. Either way, we calculate it twice, and the calculation seems
        # kinda expensive, so I'm defining it here to avoid that issue.
        dist_list = list(map(lambda y: _unknownhelper_1(state, x, maxMOI, chosenlocus, allpossiblerecrud, y), np.arange(0, allpossiblerecrud.shape[0])))

        closestrecrud = np.argmin(dist_list)
        state.mindistance[x, chosenlocus] = abs(state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[0][closestrecrud]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[1][closestrecrud]])
        state.alldistance[x, chosenlocus, 0:allpossiblerecrud.shape[0]] = dist_list
        state.allrecrf[x, chosenlocus, 0:allpossiblerecrud.shape[0]] = state.recodedf[x, maxMOI * chosenlocus + allpossiblerecrud[1]]
        state.recr0[x, chosenlocus] = maxMOI * chosenlocus + allpossiblerecrud[0][closestrecrud]
        state.recrf[x, chosenlocus] = maxMOI * chosenlocus + allpossiblerecrud[1][closestrecrud]
    else:
        old = state.recoded0[x, chosen].astype(np.int64)
        new = rand.choice(np.arange(state.frequencies_RR[0][chosenlocus])) + 1
        if is_chosen_valid:
            oldalleles = state.recoded0[x, np.intersect1d(np.arange((chosenlocus * maxMOI), (chosenlocus + 1) * maxMOI), np.where(state.hidden0[x] == HiddenAlleleType.MISSING.value)[0])]
        else:
            oldalleles = state.recodedf[x, np.intersect1d(np.arange((chosenlocus * maxMOI), (chosenlocus + 1) * maxMOI), np.where(state.hidden0[x] == HiddenAlleleType.MISSING.value)[0])]

        newallele_length = (np.mean((alleles_definitions_RR[chosenlocus]["0"][new-1], alleles_definitions_RR[chosenlocus]["1"][new-1])) + rand.normal(0, state.frequencies_RR[2][chosenlocus], 1))[0]

        repeatedold = 1 if np.sum(oldalleles == old) >= 1 else state.qq
        repeatednew = 1 if np.sum(oldalleles == new) >= 1 else state.qq

        inputVectors = list(itertools.product(np.arange(state.MOIf[x]), np.arange(state.MOI0[x])))
        allpossiblerecrud = pd.DataFrame(inputVectors)
        order = [1, 0] # setting column's order
        allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
        allpossiblerecrud.columns = [0, 1]

        if is_chosen_valid:
            tempalleles = state.alleles0[x, maxMOI * chosenlocus: maxMOI * chosenlocus + maxMOI]
            tempalleles[chosen - chosenlocus * maxMOI] = newallele_length
            temprecoded = state.recoded0[x, maxMOI * chosenlocus: maxMOI * chosenlocus + maxMOI]
            temprecoded[chosen - chosenlocus * maxMOI] = new - 1

            # This was seen in the earlier branch as well, I believe
            # this is a distance calculation taking into account that we have
            # no day0 data to work with. Either way, it was calculated twice,
            # so I have put the list here again.

            dist_list = list(map(lambda y: _unknownhelper_2(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y), np.arange(0, allpossiblerecrud.shape[0])))

            newclosestrecrud = np.argmin(dist_list)
            newmindistance = abs(tempalleles[allpossiblerecrud[0][newclosestrecrud]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[1][newclosestrecrud]])
            newalldistance = dist_list
            newallrecrf = state.recodedf[x, maxMOI * chosenlocus + allpossiblerecrud[1]]
        else:
            tempalleles = state.allelesf[x, maxMOI * chosenlocus: maxMOI * chosenlocus + maxMOI]
            tempalleles[chosen - chosenlocus * maxMOI] = newallele_length
            temprecoded = state.recodedf[x, maxMOI * chosenlocus: maxMOI * chosenlocus + maxMOI]
            temprecoded[chosen - chosenlocus * maxMOI] = new - 1

            # Same deal as twice before: Repeated calculation saved to this variable here.
            dist_list = list(map(lambda y:_unknownhelper_3(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y) , np.arange(0, allpossiblerecrud.shape[0])))

            newclosestrecrud = np.argmin(dist_list)
            newmindistance = abs(tempalleles[allpossiblerecrud[1][newclosestrecrud]] - state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[0][newclosestrecrud]])
            newalldistance = dist_list
            newallrecrf = temprecoded[allpossiblerecrud[1]]

        ## likelihoodnew
        likelihoodnew_numerator = state.dvect[np.round(newalldistance).astype(np.int64)]
        likelihoodnew_demominator = np.sum(state.frequencies_RR[1][chosenlocus][:state.frequencies_RR[0][chosenlocus]] * state.dvect[state.correction_distance_matrix[chosenlocus][newallrecrf[0].astype(np.int64)].astype(np.int64)])
        likelihoodnew = np.nanmean(likelihoodnew_numerator / likelihoodnew_demominator) * repeatednew

        ## likelihoodold
        temp = np.round(state.alldistance[x, chosenlocus])
        temp = temp[~np.isnan(temp)].astype(np.int64)
        likelihoodold_numerator = state.dvect[temp]

        ### NEED TO REVISIT
        temp_allrecrf = state.allrecrf[x, chosenlocus, :maxMOI**2]
        temp_allrecrf = temp_allrecrf[~np.isnan(temp_allrecrf)].astype(np.int64)

        temp = state.frequencies_RR[1][chosenlocus, :state.frequencies_RR[0][chosenlocus]] * state.dvect[state.correction_distance_matrix[chosenlocus][temp_allrecrf].astype(np.int64)]
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

        if z >= alpha:
            return

        state.recoded0[x, chosen] = new - 1
        if is_chosen_valid:
            state.alleles0[x, chosen] = newallele_length
        else:
            state.allelesf[x, chosen] = newallele_length
        state.mindistance[x, chosenlocus] = newmindistance
        state.alldistance[x, chosenlocus, :allpossiblerecrud.shape[0]] = newalldistance
        state.allrecrf[x, chosenlocus, :allpossiblerecrud.shape[0]] = newallrecrf
        state.recr0[x, chosenlocus] = maxMOI * (chosenlocus) + allpossiblerecrud[0][newclosestrecrud]
        state.recrf[x, chosenlocus] = maxMOI * (chosenlocus) + allpossiblerecrud[1][newclosestrecrud]

# TODO: Determine what this subroutine actually is.
# I think its for calculating a distance between microsatelite alleles. Is that right?
def _unknownhelper_1(state,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[0][y]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[1][y]])
    return value

# TODO: Determine what this subroutine actually is.
# I think its for calculating a distance between microsatelite alleles, when we don't have day 0 data. Is that right?
def _unknownhelper_2(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(tempalleles[allpossiblerecrud[0][y]] - state.allelesf[x, maxMOI * chosenlocus + allpossiblerecrud[1][y]])
    return value

# TODO: Determine what this subroutine actually is.
# Looks similar to our other evaluations, but we use day0 data. I guess this is the case where we don't have f data, but we have 0.
def _unknownhelper_3(state,tempalleles,x,maxMOI,chosenlocus,allpossiblerecrud,y):
    value = np.abs(tempalleles[allpossiblerecrud[1][y]] - state.alleles0[x, maxMOI * chosenlocus + allpossiblerecrud[0][y]])
    return value
