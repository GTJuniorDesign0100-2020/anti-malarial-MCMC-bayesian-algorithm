import pandas as pd
import numpy as np
import math

import api.recrudescence_utils as recrudescence_utils

def define_alleles(genotypedata: pd.DataFrame, locirepeats: np.ndarray, maxk: np.ndarray):
	'''
	Generate definitions of alleles (i.e. binning)
	Raw allele lengths are converted and separated into small allele classes. 
	And when it does, sometimes the data are reported as intermediate values.
	This method resolves (bins) the intermediate values.

	:param genotypedata:
			type: pandas dataframe
			description: genetic data, where first column (name 'Sample ID') has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the xth allele detected
	:param locirepeats:
			type: numpy.array
			description: a vector of length number of loci with type of locus (dinucleotide, trinucleotide, etc. repeats)
	:param maxk:
			type: numpy.array
			description: a vector of length of loci with the maximum number of alleles for each locus
	:return:
			type: list that contains dataframe
			description: list of length number of loci, each entry is a number of alleles by 2 matrix (1st column = lower bound, 2nd column = upper bound)
	'''
	
	locinames = recrudescence_utils.get_locinames(genotypedata)
	nloci = len(locinames)

	# first section
	alleles = _getAlleles(genotypedata, nloci, locirepeats, locinames)

	# second section
	compressed_alleles = _compress(alleles, nloci, maxk)

	return compressed_alleles

def _getAlleles(genotypedata: pd.DataFrame, nloci: int, locirepeats: np.ndarray, locinames: dict) -> list:
	'''
	Returns a list that contains all lower bound, upper bound, 
	and the counts of raw alleles between lower and upper bound values
	after binning the raw allele data

	:param genotypedata: The data table to retrieve allele values
	:param nloci: The number of loci
	:param locirepeats: a vector of length number of loci with type of locus
	:param locinames: The dictionary that contains all locus prefix and its columns from data table
	'''

	alleles = []
	current_column = 0
	for i in range(nloci):
		# retrieve raw alleles (each index contains every raw alleles data with the same locinames)
		# ex. all data with X313 prefix lociname in index 0
		loci_name_prefix, last_column = locinames.get(i)

		raw_alleles, current_column = recrudescence_utils.get_RawAlleles(genotypedata, current_column, last_column)

		if (max(raw_alleles) - min(raw_alleles)) < locirepeats[i]:
			# find break values(lower and upper)
			lower_break_value, upper_break_value, counts_column = _findBreakValues(raw_alleles, locirepeats, i)
			alleles = _prepareDataframes(alleles, lower_break_value, upper_break_value, counts)
		else:
			lower_break_value, upper_break_value, counts = _binRawAlleleValues(raw_alleles, locirepeats, i)
			alleles = _prepareDataframes(alleles, lower_break_value, upper_break_value, counts)

	return alleles

def _findBreakValues(raw_alleles: list, locirepeats: np.ndarray, i: int):
	'''
	Find the correct lower bound value and upper bound value given raw allele values
	'''

	lower_break_value = []
	upper_break_value = []
	counts_column = []

	lower_list = []
	upper_list = []

	for allele in raw_alleles:
		lower_list.append(allele - locirepeats[i]/2)
		upper_list.append(allele + locirepeats[i]/2)

	# search for the min from the lower_list and upper_list and add to break lists.
	lower_break_value.append(min(lower_list))
	upper_break_value.append(max(upper_list))
	counts_column.append(len(lower_list))
	return lower_break_value, upper_break_value, counts_column

def _binRawAlleleValues(raw_alleles: list, locirepeats: np.ndarray, i: int):
	'''
	Find the correct lower bound value and upper bound value given raw allele values
	This method happens when there is intermediate raw allele values to bin

	'''

	# making breaks (not sure if we need this)
	min_num = math.floor(min(raw_alleles)) - 0.5
	max_num = max(raw_alleles) + 1
	breaks = np.array([])
	while min_num < max_num:
		breaks = np.append(breaks, min_num)
		min_num += 1

	breaks_min = math.floor(min(breaks))
	breaks_max = math.floor(max(breaks))

	# allele values
	allele_values = np.array(np.round((np.array(breaks[1:]) + np.array(breaks[0:-1])) / 2))

	# historgram contains the frequency of occurrences for each breaks
	histogram = {(k+0.5): 0 for k in range(breaks_min, breaks_max)}
	for allele in raw_alleles:
		bound = math.floor(allele) + 0.5
		if allele > bound:
			histogram[bound] += 1
		else:
			histogram[bound-1] += 1

	# hist_alleles_count
	# list that contains the count for each break
	hist_alleles_count = list(histogram.values())

	# list that contains sum of 'count' from the hist_alleles_count
	# increment 'x' index of the hist_alleles_count by locirepeats[j] to select 'count'
	counts_by_offset = []
	for j in range(locirepeats[i]):
		seq = list(range(j, len(hist_alleles_count), locirepeats[i]))
		selected_count = []
		for num in seq:
			selected_count.append(hist_alleles_count[num])
		sum = 0
		for num in selected_count:
			sum += num
		counts_by_offset.append(sum)

	# select certain allele values based on the raw alleles, counts_by_offset
	seq = list(range(counts_by_offset.index(max(counts_by_offset)), len(allele_values), locirepeats[i]))
	possible_alleles = []
	for num in seq:
		possible_alleles.append(allele_values[num])

	if min(raw_alleles) <= (min(possible_alleles) - locirepeats[i]/2):
		possible_alleles = [min(possible_alleles) - locirepeats[i]] + possible_alleles

	if max(raw_alleles) > (max(possible_alleles) + locirepeats[i]/2):
		possible_alleles = possible_alleles + [max(possible_alleles) + locirepeats[i]]

	# assign clusters
	clusters = []
	for allele in raw_alleles:
		arr = np.array(possible_alleles) - allele
		arr = abs(arr).tolist()
		min_index = arr.index(min(arr))
		clusters.append(min_index)

	unique_clusters = list(dict.fromkeys(clusters))
	k = len(unique_clusters)

	# find break values(lower and upper)
	lower_break_value = []
	upper_break_value = []
	for cluster in unique_clusters:
		lower_break_value.append(possible_alleles[cluster] - locirepeats[i]/2)
		upper_break_value.append(possible_alleles[cluster] + locirepeats[i]/2)
	lower_break_value = sorted(lower_break_value)
	upper_break_value = sorted(upper_break_value)

	counts = []
	for j in range(len(lower_break_value)):
		sum = 0
		for allele in raw_alleles:
			if allele > lower_break_value[j] and allele <= upper_break_value[j]:
				sum += 1
		counts.append(sum)

	return lower_break_value, upper_break_value, counts

def _prepareDataframes(alleles: list, lower_break_value: list, upper_break_value: list, counts: list):
	'''
	Turn all lists of alleles information to pandas dataframe
	'''

	# prepare columns of lower_bound, upper_bound, and count
	allele_low = pd.DataFrame(lower_break_value)
	allele_high = pd.DataFrame(upper_break_value)
	allele_count = pd.DataFrame(counts)

	# put allele columns together to make dataframe
	allele_df = pd.concat([allele_low, allele_high, allele_count], axis=1)
	allele_df.columns = ['lower_break_value', 'upper_break_value', 'counts']
	alleles.append(allele_df)
	return alleles

def _compress(alleles: list, nloci: int, maxk: int) -> list:
	'''
	Compress alleles from the binned alleles list and returns maxk most frequent alleles form that alleles list
	
	:param alleles: the list of alleles after binning
	:param nloci: the number of loci names
	:param maxk: the number of loci with the maximum number of alleles for each locus
	'''
	alleles2 = []
	for i in range(nloci):
		sortedindex = []
		current_count = alleles[i]['counts'].tolist()

		if (len(current_count) <= maxk[i]):
			sortedindex = sorted(range(len(current_count)), key=lambda k: current_count[k], reverse=True)
		else:
			sortedindex = sorted(range(len(current_count)), key=lambda k: current_count[k], reverse=True)[:maxk[i]]

		# find new lower break and new upper break values
		new_lower_break_value = []
		new_upper_break_value = []
		for index in sortedindex:
			new_lower_break_value.append(alleles[i]['lower_break_value'][index])
			new_upper_break_value.append(alleles[i]['upper_break_value'][index])

		# put allele2 columns together to make dataframe
		alleles2_df = pd.DataFrame(columns=['0','1'])
		alleles2_df['0'] = new_lower_break_value
		alleles2_df['1'] = new_upper_break_value
		alleles2.append(alleles2_df)
	return alleles2