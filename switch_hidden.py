import random
import math
import numpy as np
import itertools

from mcmc import *

def switch_hidden(x, hidden0, hiddenf, classification, nloci, maxMOI, recoded0, frequencies_RR, alleles_definitions_RR, qq, dvect, correction_distance_matrix, alleles0, MOI0, MOIf, allelesf, mindistance, alldistance, recodedf, allrecrf, recr0, recrf, recr_repeats0, recr_repeatsf):
	z = random.uniform(0,1)
	if (np.nansum(np.concatenate((hidden0[x], hiddenf[x]))) > 0):
		if len(np.where(np.concatenate((hidden0[x], hiddenf[x])) == 1)[0]) > 1:
			chosen = np.random.choice(np.where(np.concatenate((hidden0[x], hiddenf[x])) == 1)[0]) + 1
		else:
			chosen = np.where(np.concatenate((hidden0[x], hiddenf[x])) == 1)[0][0] + 1

		if (classification[x] == 0): # REINFECTION
			if chosen <= (nloci * maxMOI):
				chosenlocus = math.ceil(chosen/maxMOI)
				old = recoded0[x][chosen - 1].astype(np.int64)
				new = np.random.choice(np.arange(frequencies_RR[0][chosenlocus - 1])) + 1
				oldalleles = recoded0[x, np.intersect1d(np.arange(((chosenlocus - 1) * maxMOI), chosenlocus * maxMOI), np.where(hidden0[x] == 1)[0])]
				newallele_length = (np.mean((alleles_definitions_RR[chosenlocus - 1]["0"][new-1], alleles_definitions_RR[chosenlocus - 1]["1"][new-1])) + np.random.normal(0, frequencies_RR[2][chosenlocus - 1], 1))[0]
				# newallele_length = np.mean([alleles_definitions_RR[chosenlocus]["0"][new], alleles_definitions_RR[chosenlocus]["1"][new]]) + np.random.normal(0, frequencies_RR[2][chosenlocus], 1)
				repeatedold = qq
				repeatednew = qq

				if sum(oldalleles == old) >= 1:
					repeatedold = 1
				if sum(oldalleles == new) >= 1:
					repeatednew = 1

				numerator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][new - 1].astype(np.int64)]) * repeatednew
				denominator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][old].astype(np.int64)]) * repeatednew
				if denominator != 0:
					alpha = numerator / denominator
				else:
					alpha = 0

				if z < alpha:
					recoded0[x][chosen - 1] = new - 1
					newallele_length = (np.mean((alleles_definitions_RR[chosenlocus - 1]["0"][new-1], alleles_definitions_RR[chosenlocus - 1]["1"][new-1])) + np.random.normal(0, frequencies_RR[2][chosenlocus - 1], 1))[0]
					alleles0[x][chosen - 1] = newallele_length

					inputVectors = list(itertools.product(np.arange(MOIf[x]), np.arange(MOI0[x])))
					allpossiblerecrud = pd.DataFrame(inputVectors)
					order = [1, 0] # setting column's order
					allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
					allpossiblerecrud.columns = [0, 1]
					closestrecrud = np.argmin(list(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0]))))
					mindistance[x][chosenlocus - 1] = abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][closestrecrud]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][closestrecrud]])
					alldistance[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = list(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0])))
					allrecrf[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1]]
					recr0[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][closestrecrud]
					recrf[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][closestrecrud]
					recr_repeats0[x][chosenlocus - 1] = sum(recoded0[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recr0[x][chosenlocus - 1])
					recr_repeatsf[x][chosenlocus - 1] = sum(recodedf[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recrf[x][chosenlocus - 1])
			else:
				chosen = chosen - nloci * maxMOI
				chosenlocus = math.ceil(chosen/maxMOI)
				old = recoded0[x][chosen - 1].astype(np.int64)
				new = np.random.choice(np.arange(frequencies_RR[0][chosenlocus - 1])) + 1
				oldalleles = recoded0[x, np.intersect1d(np.arange(((chosenlocus - 1) * maxMOI), chosenlocus * maxMOI), np.where(hidden0[x] == 1)[0])]
				# newallele_length = np.mean([alleles_definitions_RR[chosenlocus]["0"][new], alleles_definitions_RR[chosenlocus]["1"][new]]) + np.random.normal(0, frequencies_RR[2][chosenlocus], 1)
				repeatedold = qq
				repeatednew = qq

				if sum(oldalleles == old) >= 1:
					repeatedold = 1
				if sum(oldalleles == new) >= 1:
					repeatednew = 1

				numerator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][new - 1].astype(np.int64)]) * repeatednew
				denominator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][old].astype(np.int64)]) * repeatednew
				alpha = numerator / denominator

				if z < alpha:
					recodedf[x][chosen - 1] = new - 1
					newallele_length = (np.mean((alleles_definitions_RR[chosenlocus - 1]["0"][new-1], alleles_definitions_RR[chosenlocus - 1]["1"][new-1])) + np.random.normal(0, frequencies_RR[2][chosenlocus - 1], 1))[0]
					alleles0[x][chosen - 1] = newallele_length

					inputVectors = list(itertools.product(np.arange(MOIf[x]), np.arange(MOI0[x])))
					allpossiblerecrud = pd.DataFrame(inputVectors)
					order = [1, 0] # setting column's order
					allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
					allpossiblerecrud.columns = [0, 1]
					closestrecrud = np.argmin(list(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0]))))
					mindistance[x][chosenlocus - 1] = abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][closestrecrud]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][closestrecrud]])
					alldistance[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = list(map(lambda y: abs(alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0])))
					allrecrf[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1]]
					recr0[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][closestrecrud]
					recrf[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][closestrecrud]
					recr_repeats0[x][chosenlocus - 1] = sum(recoded0[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recr0[x][chosenlocus - 1])
					recr_repeatsf[x][chosenlocus - 1] = sum(recodedf[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recrf[x][chosenlocus - 1])
		else:
			if chosen <= (nloci * maxMOI):
				chosenlocus = math.ceil(chosen/maxMOI)
				old = recoded0[x][chosen - 1].astype(np.int64)
				new = np.random.choice(np.arange(frequencies_RR[0][chosenlocus - 1])) + 1
				oldalleles = recoded0[x, np.intersect1d(np.arange(((chosenlocus - 1) * maxMOI), chosenlocus * maxMOI), np.where(hidden0[x] == 1)[0])]
				newallele_length = (np.mean((alleles_definitions_RR[chosenlocus - 1]["0"][new-1], alleles_definitions_RR[chosenlocus - 1]["1"][new-1])) + np.random.normal(0, frequencies_RR[2][chosenlocus - 1], 1))[0]
				repeatedold = qq
				repeatednew = qq

				if sum(oldalleles == old) >= 1:
					repeatedold = 1
				if sum(oldalleles == new) >= 1:
					repeatednew = 1

				inputVectors = list(itertools.product(np.arange(MOIf[x]), np.arange(MOI0[x])))
				allpossiblerecrud = pd.DataFrame(inputVectors)
				order = [1, 0] # setting column's order
				allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
				allpossiblerecrud.columns = [0, 1]

				tempalleles = alleles0[x][maxMOI * (chosenlocus - 1): maxMOI * (chosenlocus - 1) + maxMOI]
				tempalleles[chosen - (chosenlocus - 1) * maxMOI - 1] = newallele_length
				temprecoded = recoded0[x][maxMOI * (chosenlocus - 1): maxMOI * (chosenlocus - 1) + maxMOI]
				temprecoded[chosen - (chosenlocus - 1) * maxMOI - 1] = new - 1

				newclosestrecrud = np.argmin(list(map(lambda y: abs(tempalleles[allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0]))))
				newmindistance = abs(tempalleles[allpossiblerecrud[0][newclosestrecrud]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][newclosestrecrud]])
				newalldistance = list(map(lambda y: abs(tempalleles[allpossiblerecrud[0][y]] - allelesf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][y]]), np.arange(0, allpossiblerecrud.shape[0])))
				newallrecrf = recodedf[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[1]]

				newrecr0 = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][newclosestrecrud]
				newrecrf = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][newclosestrecrud]
				newrecr_repeats0 = sum(temprecoded == temprecoded[allpossiblerecrud[0][newclosestrecrud]])
				newrecr_repeatsf = sum(recodedf[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recodedf[x][newrecrf])


				## likelihoodnew
				likelihoodnew_numerator = dvect[np.round(newalldistance).astype(np.int64)]
				likelihoodnew_demominator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus-1][newallrecrf[0].astype(np.int64)].astype(np.int64)])
				likelihoodnew = np.nanmean(likelihoodnew_numerator / likelihoodnew_demominator) * repeatednew


				## likelihoodold
				temp = np.round(alldistance[x][chosenlocus - 1])
				temp = temp[np.logical_not(np.isnan(temp))].astype(np.int64)
				likelihoodold_numerator = dvect[temp]

				### NEED TO REVISIT
				temp_allrecrf = []
				for i in range(maxMOI * maxMOI):
					item = allrecrf[x][chosenlocus-1][i]
					if not np.isnan(item):
						temp_allrecrf.append(item)
				temp_allrecrf = np.asarray(temp_allrecrf).astype(np.int64)

				temp = frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][temp_allrecrf].astype(np.int64)]
				likelihoodold_denominator = []
				for i in temp:
					likelihoodold_denominator.append(sum(i))

				likelihoodold_denominator = np.asarray(likelihoodold_denominator)

				likelihoodold = np.nanmean(likelihoodold_numerator / likelihoodold_denominator) * repeatedold

				if likelihoodnew == likelihoodold:
					alpha = 1
				else:
					if not likelihoodold == 0:
						alpha = likelihoodnew / likelihoodold
					else:
						alpha = 0

				if z < alpha:
					recoded0[x][chosen - 1] = new - 1
					alleles0[x][chosen - 1] = newallele_length
					mindistance[x][chosenlocus - 1] = newmindistance
					alldistance[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = newalldistance
					allrecrf[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = newallrecrf
					recr0[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][newclosestrecrud]
					recrf[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][newclosestrecrud]
					recr_repeats0[x][chosenlocus - 1] = sum(recoded0[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recoded0[x][recr0[x][chosenlocus - 1].astype(np.int64)])
					recr_repeatsf[x][chosenlocus - 1] = sum(recodedf[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recodedf[x][recrf[x][chosenlocus - 1].astype(np.int64)])
			else:
				chosen = chosen - nloci * maxMOI
				chosenlocus = math.ceil(chosen/maxMOI)
				old = recoded0[x][chosen - 1].astype(np.int64)
				new = np.random.choice(np.arange(frequencies_RR[0][chosenlocus - 1])) + 1
				newallele_length = (np.mean((alleles_definitions_RR[chosenlocus - 1]["0"][new-1], alleles_definitions_RR[chosenlocus - 1]["1"][new-1])) + np.random.normal(0, frequencies_RR[2][chosenlocus - 1], 1))[0]
				oldalleles = recodedf[x, np.intersect1d(np.arange(((chosenlocus - 1) * maxMOI), chosenlocus * maxMOI), np.where(hidden0[x] == 1)[0])]
				repeatedold = qq
				repeatednew = qq

				if sum(oldalleles == old) >= 1:
					repeatedold = 1
				if sum(oldalleles == new) >= 1:
					repeatednew = 1

				inputVectors = list(itertools.product(np.arange(MOIf[x]), np.arange(MOI0[x])))
				allpossiblerecrud = pd.DataFrame(inputVectors)
				order = [1, 0] # setting column's order
				allpossiblerecrud = allpossiblerecrud[[allpossiblerecrud.columns[i] for i in order]]
				allpossiblerecrud.columns = [0, 1]

				tempalleles = allelesf[x][maxMOI * (chosenlocus - 1): maxMOI * (chosenlocus - 1) + maxMOI]
				tempalleles[chosen - (chosenlocus - 1) * maxMOI - 1] = newallele_length
				temprecoded = recodedf[x][maxMOI * (chosenlocus - 1): maxMOI * (chosenlocus - 1) + maxMOI]
				temprecoded[chosen - (chosenlocus - 1) * maxMOI - 1] = new - 1

				newclosestrecrud = np.argmin(list(map(lambda y: abs(tempalleles[allpossiblerecrud[1][y]] - alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]]), np.arange(0, allpossiblerecrud.shape[0]))))
				newmindistance = abs(tempalleles[allpossiblerecrud[1][newclosestrecrud]] - alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][newclosestrecrud]])
				newalldistance = list(map(lambda y: abs(tempalleles[allpossiblerecrud[1][y]] - alleles0[x][maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][y]]), np.arange(0, allpossiblerecrud.shape[0])))
				newallrecrf = temprecoded[allpossiblerecrud[1]]

				newrecr0 = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][newclosestrecrud]
				newrecrf = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][newclosestrecrud]
				newrecr_repeats0 = sum(recoded0[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recoded0[x][newrecr0])
				newrecr_repeatsf = sum(temprecoded == temprecoded[allpossiblerecrud[1][newclosestrecrud]])

				## likelihoodnew
				likelihoodnew_numerator = dvect[np.round(newalldistance).astype(np.int64)]
				likelihoodnew_demominator = sum(frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus-1][newallrecrf[0].astype(np.int64)].astype(np.int64)])
				likelihoodnew = np.nanmean(likelihoodnew_numerator / likelihoodnew_demominator) * repeatednew


				## likelihoodold
				temp = np.round(alldistance[x][chosenlocus - 1])
				temp = temp[np.logical_not(np.isnan(temp))].astype(np.int64)
				likelihoodold_numerator = dvect[temp]

				### NEED TO REVISIT
				temp_allrecrf = []
				for i in range(maxMOI * maxMOI):
					item = allrecrf[x][chosenlocus-1][i]
					if not np.isnan(item):
						temp_allrecrf.append(item)
				temp_allrecrf = np.asarray(temp_allrecrf).astype(np.int64)
				temp = frequencies_RR[1][chosenlocus - 1][0:frequencies_RR[0][chosenlocus - 1]] * dvect[correction_distance_matrix[chosenlocus - 1][temp_allrecrf].astype(np.int64)]
				likelihoodold_denominator = []
				for i in temp:
					likelihoodold_denominator.append(sum(i))
				likelihoodold_denominator = np.asarray(likelihoodold_denominator)
				likelihoodold = np.nanmean(likelihoodold_numerator / likelihoodold_denominator) * repeatedold

				if likelihoodnew == likelihoodold:
					alpha = 1
				else:
					if not likelihoodold == 0:
						alpha = likelihoodnew / likelihoodold
					else:
						alpha = 0

				if z < alpha:
					recoded0[x][chosen - 1] = new - 1
					allelesf[x][chosen - 1] = newallele_length
					mindistance[x][chosenlocus - 1] = newmindistance
					alldistance[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = newalldistance
					allrecrf[x][chosenlocus - 1][0:allpossiblerecrud.shape[0]] = newallrecrf
					recr0[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[0][newclosestrecrud]
					recrf[x][chosenlocus - 1] = maxMOI * (chosenlocus - 1) + allpossiblerecrud[1][newclosestrecrud]
					recr_repeats0[x][chosenlocus - 1] = sum(recoded0[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recoded0[x][recr0[x][chosenlocus - 1].astype(np.int64)])
					recr_repeatsf[x][chosenlocus - 1] = sum(recodedf[x][maxMOI * (chosenlocus - 1) : maxMOI * chosenlocus] == recodedf[x][recrf[x][chosenlocus - 1].astype(np.int64)])