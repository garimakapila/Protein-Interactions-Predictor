# Copyright © Garima Kapila

import Bio.SeqUtils as bsu, Bio.SubsMat.MatrixInfo as mat, pandas as pd
from Bio.SeqUtils import ProtParamData
from DataParser import *
from Matrices import *


# Column names, function for various categories to compare two resdiues with

def blosum_score_categories():
	names = ['Matching', 'Similar', 'Non-Similar']
	names = ['Sum_Blosum_' + name + '_Scores' for name in names]
	return names, blosum_score_difference

def counts_categories():
	names = ['Matching', 'Non-Matching', 'Non-Matching-Gap']
	names = [name + '_Count' for name in names]
	return names, counts_difference

# default matrices for interface residue information
# for 'molecule_categories' function
matrices = [
	ProtParamData.kd,   # hydrophobicity
	ProtParamData.Flex, # flexibility 
	ProtParamData.em,   # surface probability
	ProtParamData.ja,   # surface transfer energy
	total_accessible_surface_area,
	volume_buried_residue,
	volume_side_chains,
	polarity,
	polarizability,
	solvent_accessible_surface_area,
	net_charge_index_side_chains,
]
def molecule_categories():
	names = [
		'Hydrophobicity',
		'Flexibility',
		'Emini_Surface_Fractional_Probability',
		'Janin_Interior_To_Surface_Transfer_Energy',
		'Molecular_Weight',
		'Total_Accessible_Surface_Area',
		'Volume_Buried_Residue',
		'Volume_Side_Chains',
		'Polarity',
		'Polarizability',
		'Solvent_Accessible_Surface_Area',
		'Net_Charge_Index_Side_Chains',
	]
	names = ['Difference_' + name + '_Scores' for name in names]
	return names, molecule_difference

def secondary_structure_categories():
	names = ['Helix', 'Turn', 'Sheet']
	names = ['Difference_' + name + '_Count' for name in names]
	return names, secondary_structure_difference






##################
# Helper Methods #
##################

residues_list = 'ACDEFGHIKLMNPQRSTVWY'
ambiguous_list = 'BUXZ'
blosum_dictionary = untriangularize_dictionary(mat.blosum62)

### Functions for categories to compare two residues with ###

# blosum scores of matching, non-matching similar, and non-similar residues
def blosum_score_difference(residue1, residue2):

	# gap residues
	if residue1 == '-':
		score = average_blosum_scores[(residue2)]
	elif residue2 == '-':
		score = average_blosum_scores[(residue1)]
	# non-gap residue
	else:
		score = blosum_dictionary[(residue1, residue2)]
	
	# score is positive
	if score > 0:
		# matching residues
		if residue1 == residue2:
			return [score, 0, 0]
		# similar, non-matching residues
		return [0, score, 0]
	
	# non-similar residues
	return [0, 0, score]

# counts of matching, non-matching, and non-matching gap residues
def counts_difference(residue1, residue2):
	# non-matching gap residues
	if residue1 == '-' or residue2 == '-':
		return [0, 0, 1]
	# matching residues
	if residue1 == residue2:
		return [1, 0, 0]
	# non-matching residues
	return [0, 1, 0]

# difference between values in matrices and molecular weight
def molecule_difference(residue1, residue2):
	# get difference values from matrices
	values = [difference(residue1, residue2, matrix) for matrix in matrices]
	# get difference in molecular weight
	values += [difference_molecule_weight(residue1, residue2)]
	return values

# difference between counts of helix, turn, and sheet residues
def secondary_structure_difference(residue1, residue2):
	# get secondary_structure_fraction for both residues
	second_struct1 = secondary_structure_tuple(residue1)
	second_struct2 = secondary_structure_tuple(residue2)
	# find differences between secondary_structure_fractions
	differences = [feature1 - feature2
		for feature1, feature2 in zip(second_struct1, second_struct2)]
	# square differences
	differences = [round(difference ** 2, 3) for difference in differences]
	return differences


### Helper methods for functions to compare two residues with ###

# difference given a matrix (dictionary)
def difference(residue1, residue2, matrix):
	if residue1 == '-':
		score = matrix[residue2] - matrices_averages[frozenset(matrix.items())]
	elif residue2 == '-':
		score = matrix[residue1] - matrices_averages[frozenset(matrix.items())]
	else:
		score = matrix[residue1] - matrix[residue2]
	return round(score ** 2, 3)

# difference in residue molcule weights
def difference_molecule_weight(residue1, residue2):
	weight = lambda residue: bsu.molecular_weight(residue, 'protein')
	if residue1 == '-':
		score = weight(residue2) - average_molecule_weight
	elif residue2 == '-':
		score = weight(residue1) - average_molecule_weight
	else:
		score = weight(residue1) - weight(residue2)
	return round(score ** 2, 3)

# tuple of whether residue is helix, turn, or sheet
def secondary_structure_tuple(residue):
	helix = 'VIYFWL'
	turn = 'NPGS'
	sheet = 'EMAL'
	total = float(len(residues_list))
	if residue == '-':
		return (len(helix)/total, len(turn)/total, len(sheet)/total)
	return (residue in helix, residue in turn, residue in sheet)


### Estimate of residue and gap difference = average difference ###

# for each residue, calculate average blosum score between other residues
def average_value_blosum_score_dict():
	averages = {}
	for residue in residues_list:
		averages[residue] = 0
	for residues, value in blosum_dictionary.iteritems():
		residue1, residue2 = residues
		if residue1 not in ambiguous_list:
			averages[residue1] += blosum_dictionary[(residue1, residue2)]
	for residue in residues_list:
		averages[residue] /= len(residues_list)
	return averages

def average_molecule_weight():
	weights = [bsu.molecular_weight(r, 'protein') for r in residues_list]
	return sum(weights)/len(weights)

def add_average_values(dictionary_list):
	return {
		frozenset(matrix.items()): average_value_of_dict(matrix)
		for matrix in dictionary_list
	}

average_blosum_scores = average_value_blosum_score_dict()
average_molecule_weight = average_molecule_weight()
matrices_averages = add_average_values(matrices)


