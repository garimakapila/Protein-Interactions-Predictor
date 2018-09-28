# Copyright © Garima Kapila

from DataParser import *
import time


### Compare sequences ###

"""
Return names and functions for getting information, name_prefix appended
ex. 
categories = [counts_categories(), blosum_score_categories()]
get_names_and_functions(categories, name_prefix='1')
returns:
names = ['1Matching_Count', '1Matching-Gap_Count', '1Non-Matching_Count',
	'1Non-Matching-Gap_Count', '1Sum_Blosum_Matching_Scores',
	'1Sum_Blosum_Similar_Scores', '1Sum_Blosum_Non-Similar_Scores']
and
functions = counts_differences, blosum_score_differences
"""
def get_names_and_functions(categories, name_prefix=''):
	names = flatten([name for name, function in categories])
	names = [name_prefix + '_' + name for name in names]
	functions = [function for name, function in categories]
	return names, functions

"""
Get differences between global alignment sequences
functions is a list of functions that each take in 2 residues as arguments
returns combined list of differences from those functions
"""
def get_differences(sequence1, sequence2, functions):
	# intialize differences as 0
	length = len(flatten([func('A', 'A') for func in functions]))
	results = [0] * length
	# get differences for each pair of residues in sequences
	for residue_pair in zip(sequence1, sequence2):
		# update counts in results
		comparisons = flatten([func(*residue_pair) for func in functions])
		results = add_lists(results, comparisons)

	return results



# return the sum of identical non-gap residues between seq1 and seq2
def num_identical(seq1, seq2):
	return sum(1 for r1, r2 in zip(seq1, seq2) if r1 == r2 and r1 != '-')



### Get interface residues ###

"""
increment values of a list from start_index to end_index inclusive
ex. shift_gap_indices('--A-B--CDE-', [1, 2, 3, 4, 5]) = [3, 5, 8, 9, 10]
"""
def shift_gap_indices(sequence, unzipped_indices):
	for sequence_index, residue in enumerate(sequence):
		if residue == '-':
			unzipped_indices = [
				index + 1 if index > sequence_index
				else index
				for index in unzipped_indices
			]
	return unzipped_indices


# precondition: Indices are unzipped and decremented
# Difference between protein A and B at specified indices
def get_site_differences(A, B, indices, global_alignments, functions):
	seqA, seqB = global_alignments['alnA'], global_alignments['alnB']
	residuesA, residuesB = extract_residues(seqA, seqB, indices, A)
	differences = get_differences(residuesA, residuesB, functions)
	return differences


"""
extract residues from global alignments
use sequenceA to shift indices with gaps, extract residues from sequenceA
and sequenceB at those indices
"""
def extract_residues(sequenceA, sequenceB, indices, protein):
	indices = shift_gap_indices(sequenceA, indices)
	new_indices = [i for i in indices if i < len(sequenceA)]
	if indices != new_indices:
		print '\rInvalid residue indices for %s, using a valid subset'        \
			' instead' % protein
	residuesA = extract_characters(sequenceA, new_indices)
	residuesB = extract_characters(sequenceB, new_indices)
	return residuesA, residuesB


def set_gap_dict(gap_dict):
	global gap_pattern_dict
	gap_pattern_dict = gap_dict


# indices are decremented
def gap_pattern_score(indices, global_alignments):
	seqA, seqB = global_alignments['alnA'], global_alignments['alnB']
	score = 0
	for i in indices:
		if i > 0 and i < len(seqA) - 1:
			resA_left = seqA[i - 1]
			resA = seqA[i]
			resA_right = seqA[i + 1]
			resB = seqB[i]
			if resA == '-' and (resA_left and resA_right) in gap_pattern_dict:
				score += gap_pattern_dict[(resA_left, resA_right)][resB]
	return score


"""

"""
def fragments_score():
	pass
