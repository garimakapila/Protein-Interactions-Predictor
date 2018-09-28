# Copyright © Garima Kapila

import bs4, csv, os, subprocess, sys, urllib, urllib2


"""
Fetch an organism's:
1. Fasta File from UniProt
2. Interactome from HINT
3. Interactome Interface Residue Indices from Interactome INSIDER
4. Indices of Special Sites from UniProt
5. Domain Information from The European Bioinformatics Institute (Pfam)

Database websites used:
1. UniProt: http://www.uniprot.org
2. HINT: http://hint.yulab.org
3. Interactome INSIDER: http://interactomeinsider.yulab.org
4. The European Bioinformatics Institute: https://www.ebi.ac.uk
"""


class FileFetcher():

	# input organism is of type Organism
	def __init__(self, organism):
		self.organism = organism



	# Fasta file from UniProt database
	def fetch_fasta(self, file=''):

		if file == '':
			file = self.organism.name +'.fasta'

		# print progress message
		message = 'Fetching Fasta for %s' % self.organism.info
		sys.stdout.write('\r' + message); sys.stdout.flush()

		# get page contents of url in file
		contents = self.get_fasta_page_contents()
		# no proteins exist
		if len(contents) == 0 and self.organism.strain != '':
			# try with less specific strain
			new_strain = self.organism.strain[:-1]

			# print progress message
			message = '\n> Searching Less Specific Strain: %s' % new_strain
			sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

			# get page contents of url in file
			params = self.organism.genus, self.organism.epithet
			params = params + (new_strain,)
			contents = self.get_fasta_page_contents(params)

		with open(file, 'w') as fasta:
			fasta.write(contents)
		self.format_fasta()

		return file



	# Interfaces file from Interactome INSIDER database
	def fetch_interfaces(self, file='', option='ALL'):
		
		if file == '':
			file = self.organism.name + '-Interfaces' + option.upper() + '.tsv'

		# print progress message
		message = 'Fetching Interfaces for %s' % self.organism.info
		sys.stdout.write('\r' + message); sys.stdout.flush()

		# get url to fetch interfaces from
		name = self.organism.genus[0] + '_' + self.organism.epithet.lower()
		url = 'http://interactomeinsider.yulab.org/downloads/interfaces'      \
			'{0}/{1}_interfaces{0}.txt'.format(option, name)
		urllib.urlretrieve(url, file)
		
		# find # of interactions, print progress message
		numLines = subprocess.check_output(['wc', '-l', file]).split(' ')[-2]
		numInteractions = int(numLines) - 1
		message = 'Fetched Interfaces for {}: {} interactions'                \
			.format(self.organism.info, numInteractions)
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

		return file



	# High-quality interactome binary file from HINT database
	def fetch_interactome(self, file=''):

		if file == '':
			file = self.organism.name + '-Interactome.tsv'

		# print progress message
		message = 'Fetching Interactome for %s' % self.organism.info
		sys.stdout.write('\r' + message); sys.stdout.flush()

		# fetch high-quality interactome binary file
		name = ''.join(self.organism.params)
		url = 'http://hint.yulab.org/download/%s/binary/hq/' % name
		urllib.urlretrieve(url, file)

		# find number of interactions, print progress message
		numLines = subprocess.check_output(['wc', '-l', file]).split(' ')[-2]
		numInteractions = int(numLines) - 1
		message = 'Fetched Interactome for {}: {} interactions'               \
			.format(self.organism.info, numInteractions)
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

		return file



	# Special annotated sites from UniProt database
	def fetch_special_sites(self, file='', special_sites=[]):
		
		if file == '':
			file = self.organism.name + '-Special-Sites.csv'

		# print progress message
		message = 'Fetching Special Sites for %s' % self.organism.info
		sys.stdout.write('\r' + message); sys.stdout.flush()
		
		# default special sites to gather
		if special_sites == []:
			special_sites = [
				'active site',
				'binding site',
				'calcium binding', 
				'dna binding',
				'metal binding',
				'np bind',
				'site'
			]

		# get page contents of url in file
		contents = self.get_fasta_page_contents()
		params = self.organism.params
		# no proteins exist
		if len(contents) == 0 and self.organism.strain != '':
			# try with less specific strain
			new_strain = self.organism.strain[:-1]

			# print progress message
			message = '\n> Searching Less Specific Strain: %s' % new_strain
			sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

			# get page contents of url in file
			params = self.organism.genus, self.organism.epithet
			params = params + (new_strain,)

		# fetch indices of the special sites
		url = self.get_special_sites_url(special_sites, params)
		soup = str(bs4.BeautifulSoup(urllib2.urlopen(url), 'html.parser'))
		
		with open(file, 'w') as special_sites_file:
			special_sites_file.write(soup)
		
		self.format_special_sites(file)

		message = 'Formatted Special Sites for %s' % self.organism.info
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()
		
		return file



	# Domain indices from The European Bioinformatics Institute/Pfam database
	def fetch_pfam(self, file=''):

		if file == '':
			file = self.organism.name + '-Pfam.tsv'

		# print progress message
		message = 'Fetching Domain and Family Info for %s' % self.organism.info
		sys.stdout.write('\r' + message); sys.stdout.flush()
		
		pfam_file = str(self.organism.get_taxon_ID()) + '.tsv.gz'
		url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/'              \
			'Pfam31.0/proteomes/' + pfam_file
		urllib.urlretrieve(url, pfam_file)
		subprocess.call(['gunzip', pfam_file])

		os.rename(pfam_file[:-3], file)

		# format header
		with open(file, 'r') as reader:
			data = reader.read().splitlines(True)[2:]
			total = len(data)
		with open(file, 'w') as writer:
			column_names = data[0].split('<')
			character_list = '#>'
			for i, col in enumerate(column_names):
				for c in character_list:
					col = col.replace(c, '')
					column_names[i] = col.strip(' ')
			data[0] = '\t'.join(column_names)[1:]
			writer.writelines(data)

		message = 'Fetched Domain and Family Info for %s' % self.organism.info
		message += ': %d proteins' % total
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

		return file









	###################
	#  Helper methods #
	###################
	
	def get_fasta_page_contents(self, params=()):
		url = self.get_fasta_url(params)
		soup = str(bs4.BeautifulSoup(urllib2.urlopen(url), 'html.parser'))
		return soup


	def format_fasta(self, file=''):

		if file == '':
			file = self.organism.name +'.fasta'

		# print progress message
		info = self.organism.info
		message = 'Formatting Fasta for %s' % info
		sys.stdout.write('\r' + message); sys.stdout.flush()

		# read lines of file as list
		fasta = list(csv.reader(open(file)))

		fastaList, seq = [], ''
		for entry in fasta:
			line = entry[0]
			# line is protein ID, indicated by the character '&gt' = '>'
			if '&gt' in line:
				# next protein ID found, add/clear sequence of previous protein
				fastaList.append([seq])
				seq = ''
				# format and add protein ID
				name = line.split('\t')[0]
				protID = name.split('|')[1]
				line = '>' + protID
				fastaList.append([line])
			# build sequence, keep adding lines until next protein ID is found
			else:
				seq += line
		# add last protein's sequence
		fastaList.append([seq])

		# write formatted fasta into file without header
		with open(file, 'w') as outfile:
			csv.writer(outfile).writerows(fastaList[1:])
		
		numProt = len(fastaList)/2
		message = 'Formatted Fasta for {}: {} proteins'.format(info, numProt)
		sys.stdout.write('\r' + message + '\n'); sys.stdout.flush()

		return file


	def get_fasta_url(self, params=()):

		if params == ():
			params = self.organism.params

		# adjust url if organism strain parameter is given
		is_strain = (self.organism.strain != '')
		if is_strain:
			url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes&'         \
			    'query={}%20{}%20strain%20{}&format=fasta'.format(*params)
		else:
			url = 'http://www.uniprot.org/uniprot/?query=reviewed:yes%20AND'  \
				'%20organism:%22{}%20{}{}%22&format=fasta'.format(*params)
		return url

	
	def get_special_sites_url(self, special_sites, params=()):
		
		if params == ():
			params = self.organism.params

		sites = [
			'feature(' + site.split(' ')[0] + '%20' + site.split(' ')[1] + ')' 
			if len(site.split(' ')) > 1
			else 'feature(' + site
			for site in special_sites
		]

		# adjust url if organism strain parameter is given
		is_strain = (self.organism.strain != '')
		if is_strain:
			url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes&'         \
			    'query={}%20{}%20strain%20{}&format=tab'.format(*params)
		else:
			url = 'http://www.uniprot.org/uniprot/?&query=reviewed:%20yes%20' \
				'AND%20organism:%20{}%20{}%20{}&format=tab'.format(*params)
		url += '&columns=id,{}'.format(','.join(tuple(sites)))

		return url


	# group all indices across special sites together
	def format_special_sites(self, file):
		rows = list(csv.reader(open(file), delimiter='\t'))

		for i, row in enumerate(rows[1:]):
			for j in range(1, len(row)):
				row[j] = row[j].split(' ')
				k = 0
				res = '['
				while k < len(row[j]):
					if row[j][k].isdigit():
						num1 = int(row[j][k])
						k += 1
						if row[j][k].isdigit():
							num2 = int(row[j][k])
							if len(res) > 1:
								res += ','
							if num1 == num2 and num1: res += str(num1) + ''
							else: res += str(num1) + '-' + str(num2) + ''
					k += 1
				rows[i+1][j] = res + ']'

		with open(file, 'w') as outfile:
			csv.writer(outfile).writerows(rows)


