from DataCollector.BlastP import *
from DataCollector.DataParser import *
from DataCollector.FileFetcher import FileFetcher
from DataCollector.GlobalAlignment import *
from DataCollector.Interologs import *
from DataCollector.Orthologs import *
from Organism import Organism
import shutil


class Collector():

	# input organism is of type Organism
	def __init__(self, organism1, organism2, option):
		self.organism1 = organism1
		self.organism2 = organism2
		self.names = organism1.name + '_' + organism2.name
		self.option = option.upper()
		
		# Put all results files in directory with prefix 'Data'
		self.directory = self.option + '/Data/'

		# directories to check if blast file already exists
		self.reverse_names = organism2.name + '_' + organism1.name

		# directory to move collected files to
		self.directory += organism2.name + '/' + organism1.name
		if not os.path.exists(self.directory):
			os.makedirs(self.directory)
		


	# update = only features have been updated, no sequence alignment again
	def run(self, update=False):
		
		# 1. fetch fasta, interactome, and interfaces files
		fasta1, interactome1, interfaces, sites, pfam =                       \
			self.fetch_files1(self.organism1)
		fasta2, interactome2 = self.fetch_files2(self.organism2)
		
		# 2. filter
		self.filter_relevant_proteins(fasta1, interactome1, interfaces)
		self.filter_fasta(fasta2)
		
		# 3. run blast between both organisms
		blast = self.locate_file('BlastP', reverse = True)
		if blast == None and not update:
			blast = run_blastp(fasta1, fasta2)
		elif not update:
			blast = self.copy_reverse_blast_file(blast)
			print 'Using existing blast file'
		
		if not update:
		# 4. filter blast for orthologs
			orthologs = get_orthologs(blast, fasta1, fasta2)
			
		# 5. get global alignments between orthologs
			global_alignments = get_global_alignments(orthologs, fasta1,
				fasta2, interfaces)
			
		# 6. label orthologs
			label_orthologs_global_alignments(global_alignments, orthologs)
			label_orthologs_special_sites(sites, orthologs, global_alignments)
			label_orthologs_pfam(pfam, orthologs, global_alignments)
		
		else:
			orthologs = self.locate_file('Orthologs')
			global_alignments = self.locate_file('Global_Alignments')

		# 7. join pairs of orthologs and label if they are interologs or not
		interologs = label_interologs(orthologs, interactome1, interactome2)
		
		# 8. add interface residue infomration
		interologs = add_interface_information(interologs, interfaces,
			global_alignments, sites, pfam, fasta1)
		
		# 9. delete/move certain files and return path to interologs file
		self.clean()

		

	##################
	# Helper methods #
	##################

	# delete fetched files and move data files into data subdirectory and
	# return interologs file
	def clean(self):
		delete_keywords = ['.fasta', '-Interactome', '-Interfaces', 'Pfam',
			'Special-Sites']
		move_keywords = ['.csv']
		for file in os.listdir('.'):
			if any(word in file for word in delete_keywords):
				os.remove(file)
			elif any(word in file for word in move_keywords):
				shutil.move(file, self.directory + '/' + file)

	# returns copied blast file name, replaces with reversed file
	def copy_reverse_blast_file(self, blast_file_path):
		blast_file_name = 'BlastP_' + self.names + '.csv'
		shutil.copyfile(blast_file_path, blast_file_name)
		blast = pd.read_csv(blast_file_name, sep = ',')
		
		rearranged_names = [
			'B', 'A', 'E-Value', 'Alignment_Length', 'Start_B', 'End_B',
			'Start_A', 'End_A', 'Bitscore', 'Identical_Count', 'Positive_Count',
			'Mismatch_Count', 'Gap'
		]
		blast = blast[rearranged_names]

		renamed_cols = {
			'A': 'B', 'B': 'A', 'Start_A': 'Start_B', 'Start_B': 'Start_A',
			'End_A': 'End_B', 'End_B': 'End_A'
		}
		blast = blast.rename(index = str, columns = renamed_cols)

		# filter blast proteins that are in fasta file
		fasta1 = fasta_as_dict(self.organism1.name + '.fasta')
		fasta2 = fasta_as_dict(self.organism2.name + '.fasta')
		blast = blast.loc[blast['A'].isin(fasta1.keys())]
		blast = blast.loc[blast['B'].isin(fasta2.keys())]

		blast.to_csv(blast_file_name, index = False)
		return blast_file_name

	# return path where pre-existing blast file is
	def locate_file(self, prefix, suffix = '.csv', reverse=False):
		file_name = prefix + '_'
		if reverse:
			file_name += self.reverse_names + suffix
		else:
			file_name += self.names + suffix
		for (dirpath, dirnames, filenames) in os.walk('.'):
			if file_name in filenames and                                     \
				(self.option in dirpath or dirpath == '.'):
				return dirpath + '/' + file_name
		return None

	# organism to transfer from, option is 'ALL' or 'HQ' (INSIDER database)
	def fetch_files1(self, organism):
		ff = FileFetcher(organism)
		fasta_file = ff.fetch_fasta()
		interactome_file = ff.fetch_interactome()
		interfaces_file = ff.fetch_interfaces(option = self.option)
		special_sites_file = ff.fetch_special_sites()
		pfam_file = ff.fetch_pfam()
		return [fasta_file, interactome_file, interfaces_file,
			special_sites_file, pfam_file]

	# organism to transfer to
	def fetch_files2(self, organism):
		ff = FileFetcher(organism)
		fasta_file = ff.fetch_fasta()
		interactome_file = ff.fetch_interactome()
		return [fasta_file, interactome_file]

	# filter proteins in fasta that are in interactome_file and interfaces_file
	def filter_relevant_proteins(self, fasta_file, interactome_file,
		interfaces_file):
		# print progress
		message = 'Filtering Fasta for %s' % fasta_file
		sys.stdout.write('\r' + message); sys.stdout.flush()
		# get protein dictionaries from interactome_file and interfaces_file
		prots1 = interact_prots_as_dict(interactome_file, 'Uniprot_A',        \
			'Uniprot_B')
		prots2 = interact_prots_as_dict(interfaces_file, 'P1', 'P2')
		# filter fasta with protein dictionaries and get # of filtered proteins
		self.filter_fasta(fasta_file, prots1)
		num_prots = self.filter_fasta(fasta_file, prots2)
		# print prgress
		message = 'Filtered %s: %d proteins' % (fasta_file, num_prots)
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

	"""
	filter proteins in fasta that are in prot_dict
	has_ambiguous = contain residues 'BUXZ' or not
	"""
	def filter_fasta(self, fasta_file, prot_dict={}, out_file='', 
		has_ambiguous=False):
		# return fasta file as dictionary of (protein, sequence) pairs
		fasta = fasta_as_dict(fasta_file)
		# filter proteins in fasta file that are also in prot_dict
		if prot_dict != {}:
			filtered = dict(('>' + prot, seq) for (prot, seq) in
				fasta.iteritems() if prot in prot_dict)
		else:
			filtered = dict(('>' + prot, seq) for (prot, seq) in
				fasta.iteritems() if prot in fasta)
		# default is to write into fasta file
		if out_file == '':
			out_file = fasta_file
		# write filtered fasta in out_file
		ambiguous = 'BUXZ'
		with open(out_file, 'w') as f:
			for prot, seq in filtered.iteritems():
				# check if sequences with residue 'U' should be kept
				if has_ambiguous or \
				not has_ambiguous and not any((r in seq) for r in ambiguous):
					f.write(prot + '\n' + seq + '\n')
		# return number of filtered proteins
		return len(filtered)



