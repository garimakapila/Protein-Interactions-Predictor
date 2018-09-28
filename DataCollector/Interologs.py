# Copyright © Garima Kapila

from SequenceParser import *
from GlobalAlignment import *
from DataParser import *
from ResidueComparator import *
import numpy as np, sys


"""
Compare residues at interface regions using global alignments
"""
def add_interface_information(interologs_file, interfaces_file1,
	global_alignments_file, special_sites_file, pfam_file, fasta_file1,
	out_file=''):

	filter_interologs(interologs_file, interfaces_file1)

	# read input files
	global_alignments = global_alignments_as_dict(global_alignments_file)
	interologs = pd.read_csv(interologs_file, sep=',')
	interface_indices = interface_indices_as_dict(interfaces_file1)
	columns, functions = default_interface_compare_categories()
	columns += ['Interface_Residues_Length']

	# double columns for each pair of orthologs
	columns = double_columns(columns) + ['Interface_Database']
	
	special_sites = zipped_special_sites_as_dict(special_sites_file)
	domains = domain_sites_dict(pfam_file)
	gap_pattern_dict = gap_pattern_as_dict(fasta_file1)
	set_gap_dict(gap_pattern_dict)

	total_interologs = len(interologs)
	interologs = interologs[['A1', 'A2', 'B1', 'B2']]

	values = []
	for i, (A1, A2, B1, B2) in enumerate(interologs.values):
		
		progress = 100*float(i)/total_interologs
		message = 'Adding Interface Information \t\t%d%%\t%d/%d' %            \
			(progress, i, total_interologs)
		sys.stdout.write('\r' + message); sys.stdout.flush()
			
		interface_dict = interface_indices[(A1, A2)]

		indices1 = interface_dict['indices1']
		indices2 = interface_dict['indices2']
		db = interface_dict['database']
		
		global_alns1 = global_alignments[(A1, B1)]
		global_alns2 = global_alignments[(A2, B2)]
		indices1 = parse_zipped_indices(indices1)
		indices2 = parse_zipped_indices(indices2)

		diff1 = get_site_differences(A1, B1, indices1, global_alns1, functions)
		diff2 = get_site_differences(A2, B2, indices2, global_alns2, functions)
		gap_score1 = gap_pattern_score(indices1, global_alns1)
		gap_score2 = gap_pattern_score(indices2, global_alns2)

		value_pair_1 = diff1 + [len(indices1)]
		value = value_pair_1 + diff2 + [len(indices2), db]

		overlaps1 = overlapping_interface(A1, special_sites, domains, indices1)
		overlaps2 = overlapping_interface(A1, special_sites, domains, indices2)

		value += overlaps1 + overlaps2 + [gap_score1, gap_score2]
		values.append(value)

	# print progress
	message = 'Finished Adding Interface Information \t100%\t{0}/{0}'         \
		.format(total_interologs)
	sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

	# combine dataframes
	columns += overlapping_names('Pair_1') + overlapping_names('Pair_2')
	columns += ['Interface_Gap_Score_Pair_1', 'Interface_Gap_Score_Pair_2']
	values = pd.DataFrame(values, columns = columns)
	values = values.round(3)
	interologs = pd.read_csv(interologs_file, sep=',')
	interologs = pd.concat([interologs, values], axis = 1)
	interologs = move_col_to_end(interologs, 'Database')
	interologs = move_col_to_end(interologs, 'Label')

	# write interologs
	if out_file == '':
		out_file = interologs_file
	interologs.to_csv(out_file, sep=',', index=False)

	return out_file


# Default Interface Resdues Comparision Information
def default_interface_compare_categories():
	# get counts, blosum scores, molecule, and secondary structure information
	categories = [
		counts_categories(),
		blosum_score_categories(),
		molecule_categories(),
		secondary_structure_categories()
	]
	return get_names_and_functions(categories, name_prefix='Interface')


# filter interologs that are in the interfaces file
def filter_interologs(interologs_file, interfaces_file1, out_file=''):
	interface_indices = interface_indices_as_dict(interfaces_file1)
	interologs = pd.read_csv(interologs_file, sep=',')
	A1_Values, A2_Values = interologs['A1'], interologs['A2']
	filtered = []
	for A1, A2, values in zip(A1_Values, A2_Values, interologs.values):
		if (A1, A2) in interface_indices:
			filtered.append(values)
	filtered = pd.DataFrame(filtered, columns = interologs.columns.values)
	if out_file == '':
		out_file = interologs_file
	filtered.to_csv(out_file, sep=',', index=False)


# if proteins for 1st and 2nd organism interact, label 1, else 0 for only 1st
def label_interologs(orthologs_file, interactome_fileA, interactome_fileB, 
	file_name=''):
	
	sys.stdout.write('\rLabeling Interologs'); sys.stdout.flush()

	# read input files
	orthologs = pd.read_csv(orthologs_file, sep = ',')
	interactomeA = interactome_as_dict(interactome_fileA)
	interactomeB = interactome_as_dict(interactome_fileB)

	# double columns for each pair of orthologs
	columns = list(orthologs.columns)
	columns.remove('A'), columns.remove('B')
	columns = double_columns(columns)
	columns = ['A1', 'B1', 'A2', 'B2'] + columns + ['Label']

	# get columns to iterate over
	total = len(orthologs.values)
	values = orthologs.drop(columns = ['A', 'B'], axis = 1)
	orthologs = zip(orthologs.A, orthologs.B, values.values)

	# loop through possible interolog combinations
	interologs = []
	for i, (a1, b1, values1) in enumerate(orthologs):
		for j, (a2, b2, values2) in enumerate(orthologs):
			if i > j:
				if (a1, a2) in interactomeA:
					
					# concatenate interolog pair information
					values_list1 = np.ndarray.tolist(values1)
					values_list2 = np.ndarray.tolist(values2)
					interolog = [a1, b1, a2, b2] + values_list1 + values_list2
					
					# add interaction label
					if (b1, b2) in interactomeB:
						interolog.append(1)
					else:
						interolog.append(0)

					interologs.append(interolog)

	# write interologs and return file name
	interologs = pd.DataFrame(interologs, columns = columns)
	default_file_name = orthologs_file.replace('Orthologs', 'Interologs')
	if file_name == '': file_name = default_file_name

	interologs.to_csv(file_name, sep = ',', index = False)

	sys.stdout.write('\rLabeled Interologs \n'); sys.stdout.flush()
	print '> Total:\t %d' % len(interologs)
	print '> Positive:\t %d' % len(interologs[interologs['Label'] == 1])
	print '> Unlabeled:\t %d' % len(interologs[interologs['Label'] == 0])

	return file_name


def overlapping_names(suffix):
	names = [
		'Interface_Overlapping_Special-Sites_Count_' + suffix,
		'Interface_Overlapping_Domain_Count_' + suffix
	]
	return names

# number of same interface indices with special sites and domain indices
def overlapping_interface(A, special_sites, domains, interface_indices):
	if A in special_sites:
		if special_sites[A] == None:
			special_sites_indices = []
		else:
			special_sites_indices = parse_zipped_indices(special_sites[A])
	else:
		special_sites_indices = []

	if A in domains:
		start, end = domains[A]
		domain_indices = decrement_indices(unzip_start_end(start, end))
	else:
		domain_indices = []
	
	values = [
		num_overlapping(special_sites_indices, interface_indices),
		num_overlapping(domain_indices, interface_indices)
	]
	return values


