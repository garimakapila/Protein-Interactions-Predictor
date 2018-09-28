# Copyright © Garima Kapila

import csv, os, pandas as pd, signal, subprocess, sys, time


"""
run blastp with progress messages
ex. usage:
run_blastp('Homo_Sapiens.fasta', 'Mus_Musculus.fasta')
default returns 'Blast_Homo-Sapiens_Mus-Musculus.csv' with blast results
"""
def run_blastp(fasta_file1, fasta_file2, blast_file='', output_format=10,
	args=[], names = [], max_E_Value='1e-5'):
	
	# create database using second organism's fasta file (fasta2)
	make_blastp_DB(fasta_file2)
	
	# default blast input arguments
	# see Table C1 at https://www.ncbi.nlm.nih.gov/books/NBK279684/
	if args == []:
		args = [
			'qseqid',   'sseqid', 'evalue',   'length',
			'qstart',   'qend',   'sstart',   'send',
			'bitscore', 'nident', 'positive', 'mismatch',
			'gaps'
		]
	
	# names of blast results columns
	if names == []:
		names = [
			'A', 'B', 'E-Value', 'Alignment_Length', 'Start_A', 'End_A', 
			'Start_B', 'End_B', 'Bitscore', 'Identical_Count', 'Positive_Count',
			'Mismatch_Count', 'Gap'
		]

	# if blast_file not given, use default file
	organism1 = fasta_file1.split('.')[0]
	organism2 = fasta_file2.split('.')[0]
	if blast_file == '':
		blast_file = 'BlastP_' + organism1 + '_' + organism2 + '.csv'
	
	# shell command for running blastp
	blast_results = 'blastp -query %s -db %s -out %s -evalue %s -outfmt'      \
		% (fasta_file1, fasta_file2, blast_file, max_E_Value)
	header = str(output_format) + ' ' + ' '.join(args)
	blast_command = blast_results.split(' ') + [header]

	# run blast, remove files created for database when finished
	start_blastp(blast_command, fasta_file1, blast_file)
	blast = pd.read_csv(blast_file, sep = ',', header = None, names = names)
	blast.to_csv(blast_file, index = False)
	remove_blast_database_files(fasta_file2)
	
	# return file containing blast results
	return blast_file





###################################
#  Helper methods for run_blastp  #
###################################


def make_blastp_DB(database):
	make_blast_DB_command = 'makeblastdb -in %s -dbtype prot' % (database)
	subprocess.call(make_blast_DB_command.split(' '))


def start_blastp(blast_command, fasta_file, blast_file):
	
	# begin blast process
	blast_process = subprocess.Popen(blast_command)

	# initialize progress constants
	fasta_list = list(csv.reader(open(fasta_file)))
	total_prots = len(fasta_list)/2
	current_prot_num = 0
	
	# print messages while blast process is incomplete
	while blast_process.poll() == None:
	
		# run for 5 seconds, then pause
		os.kill(blast_process.pid, signal.SIGCONT)
		time.sleep(5)
		os.kill(blast_process.pid, signal.SIGSTOP)

		# read current protein at last line of blast file
		lastLine = subprocess.check_output(['tail', '-1', blast_file])
		prot = lastLine.split(',')[0]
		
		# find current protein's position in fasta file
		for i in range(current_prot_num, total_prots):
			if fasta_list[i][0].strip('>') == prot:
				current_prot_num = i
				break
		
		# print progress message
		progress = 100*float(current_prot_num)/total_prots
		message = 'Running BLASTp \t\t\t%d%%\t%d/%d' %                        \
			(progress, current_prot_num, total_prots)
		sys.stdout.write('\r' + message); sys.stdout.flush()
		
		# resume blast process
		os.kill(blast_process.pid, signal.SIGCONT)

	message = 'BLASTp Complete \t\t100%%\t{0}/{0}\n\n'.format(total_prots)
	sys.stdout.write('\r%s\n' % message); sys.stdout.flush()


# remove files that are not needed anymore (for blastp database)
def remove_blast_database_files(database):
	blast_suffixes = ['phr', 'pin', 'psq']
	blast_files = [database + '.' + suffix for suffix in blast_suffixes]
	for file in blast_files:
		os.remove(file)

