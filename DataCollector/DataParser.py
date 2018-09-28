# Copyright © Garima Kapila

import pandas as pd, time



### Dictionaries ###

# Returns fasta file as dictionary of (protein, sequence) pairs
def fasta_as_dict(fasta_file):
	fasta_dict = {}
	with open(fasta_file) as f:
		for line in f:
			# format protein and sequence
			protein = line.strip('>\r\n')
			sequence = f.next().strip('\r\n')
			# insert (protein, sequence) pair
			fasta_dict[protein] = sequence
	return fasta_dict

# Returns dictionary of interactions from interactome_file from HINT database
def interactome_as_dict(interactome_file):
	# read interactome file
	interactome = pd.read_csv(interactome_file, sep = '\t')
	true_interaction = [1] * len(interactome.index)
	# convert interactome to dictionary will all values as true
	zipped_interactome = zip(interactome.Uniprot_A, interactome.Uniprot_B)
	zipped_interactome = zip(zipped_interactome, true_interaction)
	interactome_dict = dict(zipped_interactome)
	# reorder proteins and update dictionary
	zipped_interactome = zip(interactome.Uniprot_B, interactome.Uniprot_A)
	zipped_interactome = zip(zipped_interactome, true_interaction)
	interactome_dict.update(dict(zipped_interactome))
	return interactome_dict


"""
Get dictionary of keys protein1, protein2 with values indices1, indices2,
database that interface residues are predicted from
Has both orders of keys (protein1, protein2) and (protein2, protein1)
"""
def interface_indices_as_dict(interfaces_file):
	interfaces = pd.read_csv(interfaces_file, sep='\t')
	indices_dict = {}
	for prot1, prot2, db, indices1, indices2 in interfaces.values:
		value = {'indices1': indices1, 'indices2': indices2, 'database': db}
		indices_dict[prot1, prot2] = value
		value = {'indices1': indices2, 'indices2': indices1, 'database': db}
		indices_dict[prot2, prot1] = value
	return indices_dict


# Returns dictionary of keys proteinA, proteinB with values sequenceA, sequenceB
def global_alignments_as_dict(global_alignments_file):
	global_alignments = pd.read_csv(global_alignments_file, sep=',')
	global_alignments_dictionary = {}
	for global_alignment in global_alignments.values:
		A, B, alignment1, alignment2, score, length = global_alignment
		values = {'alnA': alignment1, 'alnB': alignment2}
		global_alignments_dictionary[A, B] = values
	return global_alignments_dictionary


"""
collect the proteins that exist in interactome file, return as dictionary
ex.
A B  corresponds to {(A, 1), (B, 1), (C, 1)}
B C
"""
def interact_prots_as_dict(interactome_file, col1, col2):
	interactome = pd.read_csv(interactome_file, sep='\t')
	# iterate over protein columns in interfaces file
	protsP1 = interactome[col1].to_dict().iteritems()
	protsP2 = interactome[col2].to_dict().iteritems()
	# collect proteins present in interfaces file and filter
	prots = dict((prot, 1) for index, prot in protsP1)
	prots.update(dict((prot, 1) for index, prot in protsP2))
	# return dictionary of proteins
	return prots


"""
get (protein, zipped residues) dictionary
zipped residues are combined occurrances
ex.
interfaces_file values:
P1 	   P1_IRES	  P2     P2_IRES
Prot1 '[1, 2, 4]' Prot2 '[3, 4, 5]'
Prot1 '[3, 5, 4]' Prot1 '[7, 8, 10]'
Prot2 '[1, 2, 4]' Prot3 '[4, 5, 6]'
returns dictionary with entries:
Prot1 '[1, 2, 4,3, 5, 4,7, 8, 10]'
Prot2 '[3, 4, 5,1, 2, 4]'
Prot3 '[4, 5, 6]
"""
def zipped_interface_residues_as_dict(interfaces_file):
	interfaces = pd.read_csv(interfaces_file, sep='\t')
	cols = (
		interfaces.P1,
		interfaces.P1_IRES,
		interfaces.P2,
		interfaces.P2_IRES
	)
	indices_dict = {}
	for protein1, indices1, protein2, indices2 in zip(*cols):
		update_indices_dict(indices_dict, protein1, indices1)
		update_indices_dict(indices_dict, protein2, indices2)
	return indices_dict


"""
get (protein, zipped residues) dictionary
zipped residues are combined across all special sites
"""
def zipped_special_sites_as_dict(special_sites_file):
	special_sites = pd.read_csv(special_sites_file, sep=',')
	prots = special_sites['Entry']
	special_sites = special_sites.drop(columns=['Entry'])
	special_sites_dict = {}
	for prot, sites in zip(prots, special_sites.values):
		combined_sites = ''
		for site in sites:
			if site != None:
				combined_sites = combine_zipped_indices(combined_sites, site)
		for i in range(2, 10):
			commas = ',' * (10 - i)
			combined_sites = combined_sites.replace(commas, ',')
		combined_sites = '[' + combined_sites[1:].strip(',')
		if combined_sites[-1] != ']':
			combined_sites += ']'
		special_sites_dict[prot] = combined_sites
	return special_sites_dict


# get (protein, domain range (start, end)) dictionary
def domain_sites_dict(pfam_file):
	pfam = pd.read_csv(pfam_file, sep='\t')
	domain_region = zip(pfam['envelope start'], pfam['envelope end'])
	return dict(zip(pfam['seq id'], domain_region))


"""
2 dimensions, tuple of 2 items as key of either order (ex. blosum matrix)
ex. previously only (A, B) = 1 was present, now (B, A) = 1 also present
"""
def untriangularize_dictionary(triangular_dictionary):
	square_dictionary = {}
	for key, value in triangular_dictionary.iteritems():
		key1, key2 = key
		square_dictionary[key1, key2] = value
		square_dictionary[key2, key1] = value
	return square_dictionary

# get average value across all keys
def average_value_of_dict(dictionary):
	counter = 0
	for key, value in dictionary.iteritems():
		counter += value
	return value/len(dictionary)


"""
Get dictionary of 2 residues and percents # of times that different residues
occur between them
ex. ('B', 'A') : {'C': 0.5, 'D': 0.2, 'Y': 0.3} 
"""
def gap_pattern_as_dict(fasta_file):
	fasta = fasta_as_dict(fasta_file)
	pattern_dict = {}
	random = int(time.time())
	mod_random = len(fasta) / 30
	for i, protein in enumerate(fasta):
		if (i + random) % mod_random == 0:
			sequence = fasta[protein]
			for j, r_left in enumerate(sequence):
				if j < len(sequence) - 2:
					r_right = sequence[j + 2]
					r_middle = sequence[j + 1]
					if (r_left, r_right) not in pattern_dict:
						pattern_dict[r_left, r_right] = {}
						pattern_dict[r_right, r_left] = {}
					if r_middle not in pattern_dict[r_left, r_right]:
						pattern_dict[r_left, r_right][r_middle] = 1
						pattern_dict[r_right, r_left][r_middle] = 1
					else:
						pattern_dict[r_right, r_left][r_middle] += 1
						pattern_dict[r_left, r_right][r_middle] += 1
	for key in pattern_dict:
		total = sum(pattern_dict[key].values())
		for subkey in pattern_dict[key]:
			pattern_dict[key][subkey] /= float(total)
			pattern_dict[key][subkey] = round(pattern_dict[key][subkey], 3)
	return pattern_dict


### Helper Functions ###

"""
adds indices to protein key in indices_dict or concatenates to current
indices value if protein key already exists in indices_dict 
"""
def update_indices_dict(indices_dict, protein, indices):
	if protein in indices_dict:
		current_indices = indices_dict[protein]
		combined = combine_zipped_indices(current_indices, indices)
		indices_dict[protein] = combined
	else:
		indices_dict[protein] = indices
	return indices_dict

"""
concatenates 2 zipped residues
ex. combine_zipped_indices('[10, 2, 3]', '[1-7,8]')
returns '[10, 2, 3,1-7,8]'
"""
def combine_zipped_indices(zipped_indices1, zipped_indices2):
	if zipped_indices1 == zipped_indices2:
		return zipped_indices1
	range1 = zipped_indices1.strip(']')
	range2 = zipped_indices2.strip('[')
	return range1 + ',' + range2




### Lists ###


"""
unzips and decrements indices
ex. [1,2 ,3-5,2-3, 8-11] = [0, 1, 2, 3, 4, 7, 8, 9, 10]
"""
def parse_zipped_indices(zipped_indices):
	indices = unzip_indices_string(zipped_indices)
	indices = decrement_indices(indices)
	return indices

"""
input is string of indices and returns a list of indices
ex.
'[1 ,2,8,7,3 ,4, 5-8, 15-17,10-19]'
becomes 
[1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
"""
def unzip_indices_string(zipped_indices_string):
	indices_string = zipped_indices_string.strip('[]').replace(' ', '')
	entries = indices_string.split(',')
	indices = []
	for entry in entries:
		try:
			index_range = entry.split('-')
			start, end = index_range
			indices += unzip_start_end(int(start), int(end))
		except ValueError:
			if entry != '':
				indices.append(int(entry))
	indices = sorted(set(indices))
	return indices

def unzip_start_end(start, end):
	return [i for i in range(start, end + 1)]

# ex. flatten([[1], [2], [3, 4]]) = [1, 2, 3, 4]
def flatten(inflated_list):
	return [item for sublist in inflated_list for item in sublist]

def num_overlapping(list1, list2):
	return len(set(list1).intersection(list2))

"""
add values at corresponding indices in two lists
ex. add_counts([1, 4, 6], [2, 3, 5]) = [3, 7, 11]
"""
def add_lists(list1, list2):
	return [v1 + v2 for v1, v2 in zip(list1, list2)]

# ex. [1, 2, 3] becomes [0, 1, 2]
def decrement_indices(indices_list):
	return [i - 1 for i in indices_list]

"""
increment values of a list from start_index to end_index inclusive
ex. increment_values([1, 2, 3, 4, 5], 1, 3) = [1, 3, 4, 5, 5]
"""
def increment_values(values_list, start_index, end_index):
	for index, value in enumerate(values_list):
		if index >= start_index and index <= end_index:
			values_list[index] += 1
	return values_list



# return string of characters at specific indices
def extract_characters(string, indices):
	return ''.join([string[i] for i in indices])



# add prefix to each column name in columns except for those in exceptions list
def add_prefix(columns, prefix, exceptions=[]):
	columns = [
		prefix + '_' + col if col not in exceptions else col
		for col in columns
	]
	return columns


def double_columns(columns):
	
	columns1 = [col + '_Pair_1' for col in columns]
	columns2 = [col + '_Pair_2' for col in columns]
	
	return columns1 + columns2


def move_col_to_end(dataframe, col_name_keyword):
	col_name = ''
	for col in dataframe.columns.values:
		if col_name_keyword in col:
			col_name = col
	col = dataframe[col_name]
	dataframe = dataframe.drop(columns=[col_name])
	dataframe[col_name] = col
	return dataframe


