import pandas as pd
import numpy as np
import math
import statistics

""" calculate frequencies of alleles
	inputs (parameters):
		- genotypedata: 
			type: pandas dataframe
			description: genetic data, where first column (name 'Sample ID') has the id of the sample, 
						 and rest of columns have the format nameoflocus_X, where X is the xth allele detected
		- alleles_definitions:
			type: list that contains dataframe
			description: list of length number of loci
						 each entry is a number of alleles by 2 matrix (1st column = lower bound, 2nd column = upper bound)

	output:
		type: list that that contains the following:
			- index[0]
				type: numpy array
				description: the length of each list of frequencies
			- index[1]
				type: numpy matrix
				description: the matrix (number of locinames by number of alleles) that contains the frequency values
			- index[2]
				type: numpy array
				description: mean SD of within allele length
"""
def calculate_frequencies3(genotypedata, alleles_definitions):

	# retrieve IDs from the Genotypedata
	ids = genotypedata.iloc[:,0].tolist()

	# retrieve the names of locus without duplicates
	# the list will contain the names with the order that was in the genotypedata (left to right)
	col_names = list(genotypedata.columns.values)[1:]
	prev_pos = None
	locinames = {}
	lociname_index = 0
	lociname_end_index = 0
	for name in col_names:
		res = name.split("_")
		if prev_pos == None:
			prev_pos = res[0]
		elif prev_pos != res[0]:
			locinames[lociname_index] = (prev_pos, lociname_end_index)
			prev_pos = res[0]
			lociname_index += 1
			lociname_end_index += 1
		else:
			lociname_end_index += 1
			if (lociname_end_index == len(col_names)-1):
				locinames[lociname_index] = (prev_pos, lociname_end_index)

	# retrieve length of ids and  length of the list of locus names
	nids = len(ids)
	nloci = len(locinames)


	frequencies = []
	variability = []

	n = 0
	for j in range(nloci):

		# retrieve raw alleles (each index contains every raw alleles data with the same locinames)
		# ex. all data with X313 prefix lociname in index 0
		loci_name_prefix, last_index = locinames.get(j)
		raw_alleles = []
		while (n <= last_index):
			raw_alleles += genotypedata.iloc[:, n+1].tolist()
			n += 1
		raw_alleles = [loci for loci in raw_alleles if str(loci) != 'nan']


		# lower = list of lower bound values
		# high = list of upper bound values
		low = alleles_definitions[j][0]
		high = alleles_definitions[j][1]

		# length of the lower bound and upper bound list
		nrows = len(alleles_definitions[j])

		sum_list = [] # needed for storing frequency values
		sd_list = [] # standard deviation
		for i in range(nrows):
			tf_table = []
			sum = 0
			for allele in raw_alleles:
				eval = allele > low[i] and allele <= high[i]
				tf_table.append(eval)
				if eval:
					sum += 1
			sum_list.append(sum)
			
			true_items = []
			for eval_i in range(len(tf_table)):
				if tf_table[eval_i]:
					true_items.append(raw_alleles[eval_i])

			if len(true_items) > 1:
				sd_list.append(statistics.stdev(true_items))
		sum_list = np.array(sum_list)
		frequencies.append(sum_list)
		
		# mean of standard deviation
		meanSD = 0
		if (len(sd_list) > 0):
			meanSD = np.mean(sd_list)

		variability.append(meanSD)
		frequencies[j] = frequencies[j] / len(raw_alleles)

	# switch freq_length and variability list to numpy array
	freq_length = np.asarray([len(frequencies[j]) for j in range(len(frequencies))])
	variability = np.asarray(variability)

	ncol = max(freq_length)

	# initialize frequency matrix with zeros
	freqmatrix = np.zeros([nloci, ncol])
	
	# fill out each box in the frequency matrix with frequency values from frequencies list
	for j in range(nloci):
		for i in range(len(frequencies[j])):
			freqmatrix[j][i] = frequencies[j][i]

	# ret = list
	# add freq_length, frequency matrix, and variability
	# ret needed to be list type to add each numpy array and matrix
	ret = []
	ret.append(freq_length)
	ret.append(freqmatrix)
	ret.append(variability)
	return ret
