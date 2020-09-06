import pandas as pd
import numpy as np
import math

""" generate definitions of alleles (i.e. binning)
	inputs (parameters):
		- genotypedata:
			type: pandas dataframe
			description: genetic data, where first column (name 'Sample ID') has the id of the sample,
						 and rest of columns have the format nameoflocus_X, where X is the xth allele detected
		- locirepeats:
			type: numpy.array
			description: a vector of length number of loci with type of locus (dinucleotide, trinucleotide, etc. repeats)
		- maxk:
			type: numpy.array
			description: a vector of length of loci with the maximum number of alleles for each locus

	output:
			type: list that contains dataframe
			description: list of length number of loci
						 each entry is a number of alleles by 2 matrix (1st column = lower bound, 2nd column = upper bound)
"""
def define_alleles(genotypedata, locirepeats, maxk):

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

	# first section
	alleles = []
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

		if (max(raw_alleles) - min(raw_alleles)) < locirepeats[j]:
			# find break values(lower and upper)
			lower_break_value = []
			upper_break_value = []
			counts_column = []

			lower_list = []
			upper_list = []

			for allele in raw_alleles:
				lower_list.append(allele - locirepeats[j]/2)
				upper_list.append(allele + locirepeats[j]/2)

			# search for the min from the lower_list and upper_list and add to break lists.
			lower_break_value.append(min(lower_list))
			upper_break_value.append(max(upper_list))
			counts_column.append(len(lower_list))

			# prepare columns of lower_bound, upper_bound, and count
			allele_low = pd.DataFrame(lower_break_value)
			allele_high = pd.DataFrame(upper_break_value)
			allele_count = pd.DataFrame(counts)

			# put allele columns together to make dataframe
			allele_df = pd.concat([allele_low, allele_high, allele_count], axis=1)
			allele_df.columns = ['lower_break_value', 'upper_break_value', 'counts']
			alleles.append(allele_df)
		else:
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
			for i in range(locirepeats[j]):
				seq = list(range(i, len(hist_alleles_count), locirepeats[j]))
				selected_count = []
				for num in seq:
					selected_count.append(hist_alleles_count[num])
				sum = 0
				for num in selected_count:
					sum += num
				counts_by_offset.append(sum)

			# select certain allele values based on the raw alleles, counts_by_offset
			seq = list(range(counts_by_offset.index(max(counts_by_offset)), len(allele_values), locirepeats[j]))
			possible_alleles = []
			for num in seq:
				possible_alleles.append(allele_values[num])

			if min(raw_alleles) <= (min(possible_alleles) - locirepeats[j]/2):
				possible_alleles = [min(possible_alleles) - locirepeats[j]] + possible_alleles

			if max(raw_alleles) > (max(possible_alleles) + locirepeats[j]/2):
				possible_alleles = possible_alleles + [max(possible_alleles) + locirepeats[j]]

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
				lower_break_value.append(possible_alleles[cluster] - locirepeats[j]/2)
				upper_break_value.append(possible_alleles[cluster] + locirepeats[j]/2)
			lower_break_value = sorted(lower_break_value)
			upper_break_value = sorted(upper_break_value)

			counts = []
			for i in range(len(lower_break_value)):
				sum = 0
				for allele in raw_alleles:
					if allele > lower_break_value[i] and allele <= upper_break_value[i]:
						sum += 1
				counts.append(sum)


			# prepare columns of lower_bound, upper_bound, and count
			allele_low = pd.DataFrame(lower_break_value)
			allele_high = pd.DataFrame(upper_break_value)
			allele_count = pd.DataFrame(counts)

			# put allele columns together to make dataframe
			allele_df = pd.concat([allele_low, allele_high, allele_count], axis=1)
			allele_df.columns = ['lower_break_value', 'upper_break_value', 'counts']
			alleles.append(allele_df)

	# second section
	# compress
	# take maxk most frequent alleles
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
