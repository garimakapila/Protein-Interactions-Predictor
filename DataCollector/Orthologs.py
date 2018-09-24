import csv, pandas as pd, sys
from DataParser import *
from SequenceParser import *
from ResidueComparator import *



# compare # residues that match or are gaps at domain sites
def label_orthologs_pfam(pfam_file, orthologs_file,
	global_alignments_file, file_name='', residue_compare_categories=()):
	sys.stdout.write('\rLabeling Domain Sites'); sys.stdout.flush()

	# read input files
	orthologs = pd.read_csv(orthologs_file, sep = ',')
	domains = domain_sites_dict(pfam_file)
	global_alignments = global_alignments_as_dict(global_alignments_file)

	# get default information categories for comparing residues
	if residue_compare_categories == ():
		residue_compare_categories = default_domain_compare_categories()
	names, functions = residue_compare_categories

	results = []
	for A, B in zip(orthologs['A'], orthologs['B']):
		global_alns = global_alignments[(A, B)]
		if A in domains:
			start, end = domains[A]
			indices = decrement_indices(unzip_start_end(start, end))
		else:
			indices = []
		info = get_site_differences(A, B, indices, global_alns, functions)
		results.append(info + [len(indices)])

	# combine global alignment values with result values
	results = pd.DataFrame(results, columns = names + ['Domain_Length'])
	results.round(3)
	orthologs = pd.concat([orthologs, results], axis = 1)

	# write global alignments into file and return file name
	if file_name == '':
		file_name = orthologs_file
	orthologs.to_csv(file_name, sep = ',', index = False)
	sys.stdout.write('\rLabeled Domain Sites \n'); sys.stdout.flush()

	return file_name



# compare and get # residues that match or are gaps at special sites
def label_orthologs_special_sites(special_sites_file, orthologs_file,
	global_alignments_file, file_name='', residue_compare_categories=()):
	sys.stdout.write('\rLabeling Special Sites'); sys.stdout.flush()

	# read input files
	orthologs = pd.read_csv(orthologs_file, sep = ',')
	global_alignments = pd.read_csv(global_alignments_file, sep = ',')
	special_sites = zipped_special_sites_as_dict(special_sites_file)
	global_alignments = global_alignments_as_dict(global_alignments_file)

	# get default information categories for comparing residues
	if residue_compare_categories == ():
		residue_compare_categories = default_special_sites_compare_categories()
	names, functions = residue_compare_categories

	results = []
	for A, B in zip(orthologs['A'], orthologs['B']):
		global_alns = global_alignments[(A, B)]
		if A in special_sites:
			if special_sites[A] == None:
				indices = []
			else:
				indices = parse_zipped_indices(special_sites[A])
		else:
			indices = []
		info = get_site_differences(A, B, indices, global_alns, functions)
		results.append(info + [len(indices)])

	# combine global alignment values with result values
	results = pd.DataFrame(results, columns = names + ['Special-Sites_Length'])
	results.round(3)
	orthologs = pd.concat([orthologs, results], axis = 1)

	# write global alignments into file and return file name
	if file_name == '':
		file_name = orthologs_file
	orthologs.to_csv(file_name, sep = ',', index = False)
	sys.stdout.write('\rLabeled Special Sites \n'); sys.stdout.flush()

	return file_name


"""
compare and get # residues that match or are gaps in global alignment
precondition: each row in global_alignments and orthologs_file correspond
to the same set of proteins
"""
def label_orthologs_global_alignments(global_alignments_file,
	orthologs_file, file_name='', residue_compare_categories=()):
	
	sys.stdout.write('\rLabeling Global Alignments'); sys.stdout.flush()

	# read input files
	global_alignments = pd.read_csv(global_alignments_file, sep = ',')
	columns = list(global_alignments.columns.values)
	exceptions = ['A', 'B', 'Fasta_Length_A', 'Fasta_Length_B']
	columns = add_prefix(columns, 'Global', exceptions=exceptions)
	global_alignments.columns = columns

	orthologs = pd.read_csv(orthologs_file, sep = ',')
	columns = list(orthologs.columns.values)
	columns = add_prefix(columns, 'Blast', exceptions=exceptions)
	orthologs.columns = columns

	# get default information categories for comparing residues
	if residue_compare_categories == ():
		residue_compare_categories = default_global_compare_categories()
	names, functions = residue_compare_categories

	results = []
	for row in global_alignments.values:
		protA, protB, alignment1, alignment2, score, length = row
		info = get_differences(alignment1, alignment2, functions)
		results.append(info)
	results = pd.DataFrame(results, columns = names)
	results.round(3)

	# remove global alignment sequences
	global_alignments.drop('Global_Alignment1', axis = 1, inplace = True)
	global_alignments.drop('Global_Alignment2', axis = 1, inplace = True)

	# combine global alignment values with result values
	global_alignments = pd.concat([global_alignments, results], axis = 1)
	global_alignments.drop(columns = ['A', 'B'], axis = 1, inplace = True)
	orthologs = pd.concat([orthologs, global_alignments], axis = 1)

	# write global alignments into file and return file name
	if file_name == '':
		file_name = orthologs_file
	orthologs.to_csv(file_name, sep = ',', index = False)
	sys.stdout.write('\rLabeled Global Alignments \n'); sys.stdout.flush()

	return file_name


"""
input protType 'A' or 'B' in blast file, the blast file, and corresponding
fasta file of protType protein
ex. calculate_coverage('A', blast, 'HomoSapiens.fasta') where 'A' represents
Homo Sapiens protein
"""
def get_orthologs(blast_file, fasta_file1, fasta_file2, orthologs_file='',
	threshold=0.5):
	
	# print progress message
	sys.stdout.write('\rFiltering Orthologs'); sys.stdout.flush()

	# read blast and drop duplicate entries
	blast = pd.read_csv(blast_file, sep = ',')
	blast = blast.drop_duplicates(['A', 'B'])
	
	# calculate sequence coverages in blast alignment sequence
	orthologs = calculate_coverage('A', blast, fasta_file1, threshold)
	orthologs = calculate_coverage('B', orthologs, fasta_file2, threshold)

	# if orthologs_file not specified, use default file to write results in
	if orthologs_file == '':
		name = fasta_file1.split('.')[0] + '_' + fasta_file2.split('.')[0]
		orthologs_file = 'Orthologs_' + name + '.csv'
	
	# print progress message
	numProts = len(orthologs.index)
	message = 'Filtered Orthologs: %d pairs' % numProts
	sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

	# write and return orthologs_file
	orthologs.to_csv(orthologs_file, sep = ',', index = False)
	return orthologs_file



### Helper functions ###

# Default Domain Information: counts and blosum scores
def default_domain_compare_categories():
	categories = [
		counts_categories(),
		blosum_score_categories()
	]
	return get_names_and_functions(categories, name_prefix='Domain')


# Default Special Sites Information: counts and blosum scores
def default_special_sites_compare_categories():
	categories = [
		counts_categories(),
		blosum_score_categories()
	]
	return get_names_and_functions(categories, name_prefix='Special-Sites')


# Default Global Alignment Information: counts and blosum scores
def default_global_compare_categories():
	categories = [
		counts_categories(),
		blosum_score_categories()
	]
	return get_names_and_functions(categories, name_prefix='Global')


"""
input protType 'A' or 'B' in blast file, the blast file, and corresponding
fasta file of protType protein
ex. calculate_coverage('A', blast, 'HomoSapiens.fasta') where 'A' represents
Homo Sapiens protein
"""
def calculate_coverage(protType, blast, fasta_file, threshold):
	# turn 'SettingWithCopyWarning' warning off
	blast.is_copy = False
	# add column of protein sequence lengths from fasta file
	seq_length_col = 'Fasta_Length_' + protType
	fasta = fasta_as_dict(fasta_file)
	blast[seq_length_col] = blast[protType].apply(lambda prot: len(fasta[prot]))
	# columns from blast file needed to calculate coverage
	cols = [col + protType for col in ['Start_', 'End_']]
	cols.append(seq_length_col)
	# add coverage column
	coverage_col = 'Coverage_' + protType
	blast[coverage_col] = blast[cols].apply(lambda x: coverage(*x), axis = 1)
	# filter coverages by threshold and return results
	blast = blast[blast[coverage_col] >= threshold]
	# return filtered blast
	return blast

"""
helper function for calculating coverage
coverage of a protein: number residues in blast alignment/number total residues
"""
def coverage(blast_start, blast_end, fasta_seq_len):
	blastLength = blast_end - blast_start
	coverage = float(blastLength) / fasta_seq_len
	return round(coverage, 3)


