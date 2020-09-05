from random import uniform
import numpy as np
import pandas as pd
import itertools

#from calculate_frequencies import *
#from define_alleles import *
#from findposteriorfrequencies import *
#from Import_Microsattelite_Data import *
#from main import *
#from mcmc import *
#from recode_alleles import *
#from run_all_arms import *

def switch_hidden(x):
    z = random.uniform((1,))

    if sum(hidden0[x], hiddenf[x]) > 0: #if hidden alleles exist TODO: replicate na.rm
        if len(np.where(x == 1, np.concatonate(hidden0[x], hiddenf[x]))) > 1:
            chosen = np.random.choice(np.where(x == 1, np.concatonate(hidden0[x], hiddenf[x])), 1, False)
        else:
            chosen = np.where(x == 1, np.concatonate(hidden0[x], hiddenf[x]))

    if classification[x] == 0: #reinfection
        if chosen <= nloci * maxMOI: #day0 hidden allele
            chosenlocus = math.ceil(chosen / maxMOI)
            old = recoded0[x][chosen]
            new = np.random.choice(np.arange(0,frequencies_RR[0][chosenLocus]), 1, False)

            oldalleles = recoded0[x, np.arange((chosenlocus - 1) * maxMOI + 1, chosenlocus * maxMOI).intersect(np.where(x == 1, hidden0[x]))] #TODO: Double check this line, this was a rough one
            repeatedold = qq
            repeatednew = qq
            if sum(oldalleles == old) >= 1: #if old allele is a repeat, don't penalize with missing probability
                repeatedold = 1
            if sum(oldalleles == new) >= 1: #if new allele is a repeat, don't penalize with missing probability
                repeatednew = 1
            alpha = (sum(frequencies_RR[2][chosenlocus][np.arange(0, frequencies_RR[1][chosenlocus])] * dvect[correction_distace_matrix[chosenlocus][:,new] + 1])
             * repeatednew) / (sum(frequencies_RR[2][chosenlocus][np.arange(0, frequencies_RR[1][chosenlocus])] * dvect[correction_distace_matrix[chosenlocus][:,old] + 1])
             * repeatedold) #TODO: double check this line as well
            if z < alpha:
                recoded0[x][chosen] = new
                newallele_length = np.mean(alleles_definitions_RR[chosenlocus][new]) + np.random.normal(0, frequencies_RR[3][chosenlocus], 1)
                alleles0[x][chosen] = newallelelength

                #TODO: Below code had if block commented out. Double check scoping is correct. Also check closest recrud and alldistance[][][]
                inputVectors = list(itertools.product(np.arange(MOI0[x], np.arange(MOIf[x]))))
                allpossiblerecrud = pd.DateFrame(inputVectors)
                closestrecrud = np.min(np.where(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][1]]), np.arange(0, pd.shape(allpossiblerecrud)[0]))))
                mindistance[x][chosenlocus] = abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[closetrecrud][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[closetrecrud][1]])
                alldistance[x][chosenlocus][: pd.shape(allpossiblerecrud)[0]] =  map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][1]]) , np.arange(0, pd.shape(allpossiblerecrud)[0]))
                allrecrf[x][chosenlocus][: pd.shape(allpossiblerecrud)[0]] = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[:,1]]
                recr0[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[closestrecrud][0]
                recrf[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[closestrecrud][1]
                recr_repeats0[x][chosenlocus] = sum(recoded0[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus - 1) + 1)] == recoded0[x][recr0[x][chosenlocus]])
                recr_repeatsf[x][chosenlocus] = sum(recodedf[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus - 1) + 1)] == recodedf[x][recrf[x][chosenlocus]])

        #TODO: Refactor below into function
        #TODO: Check with Matt to check what to do with commented blocks in r code

        else: #day f hidden allele
            chosen = chosen - nloci * maxMOI
            chosenlocus = math.ceil(chosen / maxMOI)
            old = recodedf[x][chosen]
            new = np.random.choice(np.arange(0,frequencies_RR[0][chosenLocus]), 1, False)
            oldalleles = recodedf[x, np.arange((chosenlocus - 1) * maxMOI + 1, chosenlocus * maxMOI).intersect(np.where(x == 1, hiddenf[x]))]
            repeatedold = qq
            repeatednew = qq
            #TODO: Logic on if blocks seems very strange: confer with team and or Matt
            if sum(oldalleles == old) >= 1: #if old allele is a repeat, don't penalize with missing probability
                repeatedold = 1
            if sum(oldalleles == new) >= 1: #if new allele is a repeat, don't penalize with missing probability
                repeatednew = 1
            alpha = (sum(frequencies_RR[2][chosenlocus][np.arange(0, frequencies_RR[1][chosenlocus])] * dvect[correction_distace_matrix[chosenlocus][:,new] + 1])
             * repeatednew) / (sum(frequencies_RR[2][chosenlocus][np.arange(0, frequencies_RR[1][chosenlocus])] * dvect[correction_distace_matrix[chosenlocus][:,old] + 1])
             * repeatedold)
            if z < alpha: #switch made
                recodedf[x][chosen] = new
                newallele_length = np.mean(alleles_definitions_RR[chosenlocus][new]) + np.random.normal(0, frequencies_RR[3][chosenlocus], 1)
                allelesf[x][chosen] = newallelelength

                inputVectors = list(itertools.product(np.arange(MOI0[x], np.arange(MOIf[x]))))
                allpossiblerecrud = pd.DateFrame(inputVectors)
                closetrecrud = np.min(np.where(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][1]]), np.arange(0, pd.shape(allpossiblerecrud)[0]))))
                mindistance[x][chosenlocus] = abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[closetrecrud][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[closetrecrud][1]])
                alldistance[x][chosenlocus][: pd.shape(allpossiblerecrud)[0]] =  map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][1]]) , np.arange(0, pd.shape(allpossiblerecrud)[0]))
                allrecrf[x][chosenlocus][: pd.shape(allpossiblerecrud)[0]] = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[:,1]]
                recr0[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[closestrecrud][0]
                recrf[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[closestrecrud][1]
                recr_repeats0[x][chosenlocus] = sum(recoded0[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus - 1) + 1)] == recoded0[x][recr0[x][chosenlocus]])
                recr_repeatsf[x][chosenlocus] = sum(recodedf[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus - 1) + 1)] == recodedf[x][recrf[x][chosenlocus]])


    else: #recrudescence
        if chosen <= nloci * maxMOI: #day 0 hidden allele
            chosenlocus = math.ceil(chosen / maxMOI)
            old = recoded0[x][chosen]
            new = np.random.choice(np.arange(0,frequencies_RR[0][chosenLocus]), 1, False)
            newallele_length = np.mean(alleles_definitions_RR[chosenlocus][new]) + np.random.normal(0, frequencies_RR[3][chosenlocus], 1)
            oldalleles = recoded0[x, np.arange((chosenlocus - 1) * maxMOI + 1, chosenlocus * maxMOI).intersect(np.where(x == 1, hidden0[x]))]
            repeatedold = qq
            repeatednew = qq
            if sum(oldalleles == old) >= 1: #if old allele is a repeat, don't penalize with missing probability
                repeatedold = 1
            if sum(oldalleles == new) >= 1: #if new allele is a repeat, don't penalize with missing probability
                repeatednew = 1
            inputVectors = list(itertools.product(np.arange(MOI0[x], np.arange(MOIf[x]))))
            allpossiblerecrud = pd.DateFrame(inputVectors)
            tempalleles = alleles0[x][maxMOI * (chosenlocus - 1) + 1 : maxMOI]
            tempalleles[chosen - (chosenlocus - 1) * maxMOI] = newallelelength
            temprecoded = recoded0[x][maxMOI * (chosenlocus - 1) + 1 : maxMOI]
            temprecoded[chosen - (chosenlocus - 1) * maxMOI] = new

            newclosestrecrud = np.min(np.where(map(lambda y: abs(tempalleles[allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][1]]), np.arange(0, pd.shape(allpossiblerecrud)[0]))))
            newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud][0]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][1]])
            newalldistance = map(lambda y: abs(tempalleles[allpossiblerecrud[y][0]] - allelesf[x][maxMOI * (chosenLocus - 1) + allpossiblerecrud[y][1]]) , np.arange(0, pd.shape(allpossiblerecrud)[0]))
            newallrecrf = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[:,1]]

            #calculate new multiple-comparisons coefficient
            newrecr0 = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][0]
            newrecrf = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][1]
            newrecr_repeats0 = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud][0]]) #TODO: figure out na.rm
            newrecr_repeatsf = sum(recodedf[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus))] == recodedf[x][newrecrf])

            #TODO: Find python equivalent of na.rm
            likelihoodnew = mean(dvect[round(newalldistance) + 1] / map(lambda z: sum(frequencies_RR[1][chosenlocus][:frequencies_RR[0][chosenlocus]] * dvect[correction_distace_matrix[chosenlocus][:,newallrecrf[z]] + 1]) , np.arange(0, len(newallrecrf)))) * repeatednew
            likelihoodold = mean(dvect[round(alldistance[x][chosenlocus]) + 1] / map(lambda z: sum(frequencies_RR[1][chosenlocus][:frequencies_RR[0][chosenlocus]] * dvect[correction_distace_matrix[chosenlocus][:,allrecrf[x][chosenlocus][z]] + 1], np.arange(0, maxMOI * maxMOI)))) * repeatedold

            #TODO: Port debugging code if wanted

            if likelihoodnew == likelihoodold:
                # if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
                alpha = 1
            else:
                alpha = likelihoodnew / likelihoodold

            if z < alpha:
                #TODO: Refactor below into function
                recoded0[x][choosen] = new
                alleles0[x][choosen] = newallelelength
                mindistance[x][chosenlocus] = newmindistance
                alldistance[x][chosenlocus][:pd.shape(allpossiblerecrud)[0]] = newalldistance
                allrecrf[x][chosenlocus][:pd.shape(allpossiblerecrud)[0]] = newallrecrf
                recr0[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][0]
                recrf[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][1]
                recr_repeats0[x][chosenlocus] = sum(recoded0[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * chosenlocus)] == recoded0[x][recr0[x][chosenlocus]])
                recr_repeatsf[x][chosenlocus] = sum(recodedf[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * chosenlocus)] == recodedf[x][recr0[x][chosenlocus]])

        else: #day f hidden allele
            chosen = chosen - nloci * maxMOI
            chosenlocus = math.ceil(chosen / maxMOI)
            old = recodedf[x][chosen]
            new = np.random.choice(np.arange(0,frequencies_RR[0][chosenLocus]), 1, False)
            newallele_length = np.mean(alleles_definitions_RR[chosenlocus][new]) + np.random.normal(0, frequencies_RR[3][chosenlocus], 1)
            oldalleles = recodedf[x, np.arange((chosenlocus - 1) * maxMOI + 1, chosenlocus * maxMOI).intersect(np.where(x == 1, hidden0[x]))]
            repeatedold = qq
            repeatednew = qq
            if sum(oldalleles == old) >= 1: #if old allele is a repeat, don't penalize with missing probability
                repeatedold = 1
            if sum(oldalleles == new) >= 1: #if new allele is a repeat, don't penalize with missing probability
                repeatednew = 1
            inputVectors = list(itertools.product(np.arange(MOI0[x], np.arange(MOIf[x]))))
            allpossiblerecrud = pd.DateFrame(inputVectors)
            tempalleles = allelesf[x][maxMOI * (chosenlocus - 1) + 1 : maxMOI]
            tempalleles[chosen - (chosenlocus - 1) * maxMOI] = newallelelength
            temprecoded = recodedf[x][maxMOI * (chosenlocus - 1) + 1 : maxMOI]
            temprecoded[chosen - (chosenlocus - 1) * maxMOI] = new

            newclosestrecrud = np.min(np.where(map(lambda y: abs(tempalleles[allpossiblerecrud[y][1]] - alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[y][0]]), np.arange(0, pd.shape(allpossiblerecrud)[0]))))
            newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud][1]] - alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][0]])
            newalldistance = map(lambda y: abs(tempalleles[allpossiblerecrud[y][1]] - alleles0[x][maxMOI * (chosenLocus - 1) + allpossiblerecrud[y][0]]) , np.arange(0, pd.shape(allpossiblerecrud)[0]))
            newallrecrf = temprecoded[allpossiblerecrud][:,1]

            #calculate new multiple-comparisons coefficient
            newrecr0 = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][0]
            newrecrf = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][1]
            newrecr_repeats0 = sum(recoded0[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * (chosenlocus))] == recoded0[x][newrecrf])
            newrecr_repeatsf = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud][0]]) #TODO: check na.rm

            #TODO: Find python equivalent of na.rm
            likelihoodnew = mean(dvect[round(newalldistance) + 1] / map(lambda z: sum(frequencies_RR[1][chosenlocus][:frequencies_RR[0][chosenlocus]] * dvect[correction_distace_matrix[chosenlocus][:,newallrecrf[z]] + 1]) , np.arange(0, len(newallrecrf)))) * repeatednew
            likelihoodold = mean(dvect[round(alldistance[x][chosenlocus]) + 1] / map(lambda z: sum(frequencies_RR[1][chosenlocus][:frequencies_RR[0][chosenlocus]] * dvect[correction_distace_matrix[chosenlocus][:,allrecrf[x][chosenlocus][z]] + 1], np.arange(0, maxMOI * maxMOI)))) * repeatedold

            #TODO: Add debugging if deemed necessary

            if likelihoodnew == likelihoodold:
                # if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
                alpha = 1
            else:
                alpha = likelihoodnew / likelihoodold

            if z > alpha: #switch made
                recodedf[x][choosen] = new
                allelesf[x][choosen] = newallelelength
                mindistance[x][chosenlocus] = newmindistance
                alldistance[x][chosenlocus][:pd.shape(allpossiblerecrud)[0]] = newalldistance
                allrecrf[x][chosenlocus][:pd.shape(allpossiblerecrud)[0]] = newallrecrf
                recr0[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][0]
                recrf[x][chosenlocus] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[newclosestrecrud][1]
                recr_repeats0[x][chosenlocus] = sum(recoded0[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * chosenlocus)] == recoded0[x][recr0[x][chosenlocus]])
                recr_repeatsf[x][chosenlocus] = sum(recodedf[x][(maxMOI * (chosenlocus - 1) + 1) : (maxMOI * chosenlocus)] == recodedf[x][recr0[x][chosenlocus]])