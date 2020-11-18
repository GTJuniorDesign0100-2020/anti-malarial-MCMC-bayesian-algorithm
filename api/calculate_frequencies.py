import pandas as pd
import numpy as np
import math
import statistics

import api.recrudescence_utils as recrudescence_utils

def calculate_frequencies3(genotypedata, alleles_definitions):
	'''
	Calculate frequencies of alleles
	Using the input data table and alleles_definition 
	which contains the lower and upper break values of the allele data,
	this method calculates the amount (frequency) and the variability of raw allele data 
	in each lower and upper break values.

	Returns a list that that contains the following:
		- index[0]
			type: numpy array
			description: the length of each list of frequencies
		- index[1]
			type: numpy matrix
			description: the matrix (number of locinames by number of alleles) that contains the frequency values
		- index[2]
			type: numpy array
			description: mean SD of within allele length

	:param genotypedata
		type: pandas dataframe
		description: genetic data, where first column (name 'Sample ID') has the id of the sample,
			and rest of columns have the format nameoflocus_X, where X is the xth allele detected
	:param alleles_definitions
		type: list that contains dataframe
		description: list of length number of loci
			each entry is a number of alleles by 2 matrix (1st column = lower bound, 2nd column = upper bound)
	'''

	locinames = recrudescence_utils.get_locinames(genotypedata)
	nloci = len(locinames)

	frequencies = []
	variability = []

	n = 0
	for j in range(nloci):
		# retrieve raw alleles (each index contains every raw alleles data with the same locinames)
		# ex. all data with X313 prefix lociname in index 0
		loci_name_prefix, last_index = locinames.get(j)
		raw_alleles, n = recrudescence_utils.get_RawAlleles(genotypedata, n, last_index)

		# lower = list of lower bound values
		# high = list of upper bound values
		low = alleles_definitions[j]["0"]
		high = alleles_definitions[j]["1"]

		# length of the lower bound and upper bound list
		nrows = len(alleles_definitions[j])

		sum_list, meanSD = _get_sumList(nrows, raw_alleles, low, high)

		frequencies.append(sum_list)
		variability.append(meanSD)
		frequencies[j] = frequencies[j] / len(raw_alleles)

	# switch freq_length and variability list to numpy array
	freq_length = np.asarray([len(frequencies[j]) for j in range(len(frequencies))])
	variability = np.asarray(variability)

	ncol = max(freq_length)

	# create matrix with frequency values from frequencies list
	freqmatrix = _create_frequencyMatrix(nloci, ncol, frequencies)

	# final result
	ret = _pack_result(freq_length, freqmatrix, variability)
	return ret

def _get_sumList(nrows: int, raw_alleles: list, low: pd.core.series.Series, high: pd.core.series.Series):
	'''
	Returns a numpy array of the number of allele values that is between lower and upper bound values
	Also returns a mean of standard deviation of that numpy array.

	:param nrows: The number of rows of the alleles_definition
	:param raw_alleles: The allele values retrieved from the input file
	:param low: The list of the lower bound values from the alleles_definition
	:param high: The list of the upper bound values from the alleles_definition
	'''
	sum_list = [] # needed for storing frequency values
	sd_list = [] # standard deviation
	for i in range(nrows):
		tf_table = []
		result_sum = 0
		for allele in raw_alleles:
			eval = allele > low[i] and allele <= high[i]
			tf_table.append(eval)
			if eval:
				result_sum += 1
		sum_list.append(result_sum)

		true_items = []
		for eval_i in range(len(tf_table)):
			if tf_table[eval_i]:
				true_items.append(raw_alleles[eval_i])

		if len(true_items) > 1:
			sd_list.append(statistics.stdev(true_items))
	sum_list = np.array(sum_list)

	# mean of standard deviation
	meanSD = 0
	if (len(sd_list) > 0):
		meanSD = np.mean(sd_list)

	return sum_list, meanSD

def _create_frequencyMatrix(nloci: int, ncol: np.int64, frequencies: list):
	'''
	Turn 1D frequencies list into 2D numpy matrix

	:param nloci: The number of rows
	:param ncol: The number of columns
	:param frequencies: The 1D list that contains the frequency values
	'''

	# initialize frequency matrix with zeros
	freqmatrix = np.zeros([nloci, ncol])

	# fill out each box in the frequency matrix with frequency values from frequencies list
	for j in range(nloci):
		for i in range(len(frequencies[j])):
			freqmatrix[j][i] = frequencies[j][i]
	return freqmatrix

def _pack_result(freq_length: np.ndarray, freqmatrix: np.ndarray, variability: np.ndarray):
	'''
	Returns a list that contains all frequency-related matrices and arrays
	'''

	ret = []
	ret.append(freq_length)
	ret.append(freqmatrix)
	ret.append(variability)
	return ret