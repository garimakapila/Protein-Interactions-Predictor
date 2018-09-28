# Copyright © Garima Kapila

import Bio.pairwise2 as pairwise, Bio.SubsMat.MatrixInfo as mat, datetime, sys
from SequenceParser import *
from DataParser import *


def get_global_alignments(orthologs_file, fasta_file1, fasta_file2, 
	interfaces_file1='', file_name=''):

	date_time = datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")
	message = '\n\nStarting Global Alignments, current time: %s' % date_time
	sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

	# read input files
	fasta1 = fasta_as_dict(fasta_file1)
	fasta2 = fasta_as_dict(fasta_file2)
	orthologs = pd.read_csv(orthologs_file, sep = ',')
	total_orthologs = len(orthologs)

	# also get best global alignments for interface residues
	if interfaces_file1 != '':
		indices = zipped_interface_residues_as_dict(interfaces_file1)
	
	global_alignments = []
	
	for i, pair in enumerate(zip(orthologs.A, orthologs.B)):

		# print progress message
		progress = 100*float(i)/total_orthologs
		message = 'Getting Global Alignments \t\t%d%%\t%d/%d' %               \
			(progress, i, total_orthologs)
		sys.stdout.write('\r' + message); sys.stdout.flush()
		
		# get protein sequences
		protein_A, protein_B = pair
		sequences = fasta1[protein_A], fasta2[protein_B]
		
		# global alignment parameters if interface residues are considered
		if interfaces_file1 == '':
			args = sequences + (indices[protein_A],)
		else:
			args = sequences
		
		# get global alignment
		result = (protein_A,protein_B,) + global_alignment(*args)
		global_alignments.append(result)

	# add global alignment column names
	cols = ['A', 'B', 'Alignment1', 'Alignment2', 'Score', 'Length']
	global_alignments = pd.DataFrame(global_alignments, columns = cols)
	
	# default file name
	if file_name == '':
		organism1 = fasta_file1.split('.')[0]
		organism2 = fasta_file2.split('.')[0]
		name = organism1 + '_' + organism2
		file_name = 'Global_Alignments_' + name + '.csv'

	# print progress
	message = 'Finished Getting Global Alignments \t100%\t{0}/{0}\n\n'        \
		.format(total_orthologs)
	sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

	# write global alignments into file and return file name
	global_alignments.to_csv(file_name, sep = ',', index = False)
	return file_name


"""
Helper method for get_global_alignments
Get best alignment by number identical global, then by number identical
interface residues if zipped indices are given
"""
def global_alignment(seq1, seq2, zipped_residues='[-1]'):

	# get (optional) interface residues
	if zipped_residues != '[-1]':
		residue_indices = unzip_indices_string(zipped_residues)
		residue_indices = decrement_indices(residue_indices)
	else:
		residue_indices = []

	# assign gap penalties, more information at:
	# http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
	gap_open, gap_extend  = -10, -0.5
	global_alignment_args = (seq1, seq2, mat.blosum62, gap_open, gap_extend)
	global_alignments     = pairwise.align.globalds(*global_alignment_args)
	
	# get best alignment by number identical global, number identical interface
	best_alignment = ()
	max_global_identical, max_interface_identical = -1, -1
	
	for alignment in global_alignments:
		
		alignment1, alignment2, score, start, end = alignment
		
		# sum identical non-gap residues
		global_identical = num_identical(alignment1, alignment2)

		# sum identical non-gap residues from only the interface residues 
		if residue_indices != []:
			temp_indices = shift_gap_indices(alignment1, residue_indices)
			interface_seq1 = extract_characters(alignment1, temp_indices)
			interface_seq2 = extract_characters(alignment2, temp_indices)
			interface_identical = num_identical(interface_seq1, interface_seq2)
		else:
			interface_identical = 0

		# update maximum indentical values and best alignment
		if global_identical > max_global_identical or                         \
			(global_identical == max_global_identical and 
			interface_identical > max_interface_identical):
			
			max_global_identical = global_identical
			max_interface_identical = interface_identical
			best_alignment = alignment
	
	# return best alignment results
	alignment1, alignment2, score, start, end = best_alignment
	alignment_length = end - start

	return alignment1, alignment2, score, alignment_length

