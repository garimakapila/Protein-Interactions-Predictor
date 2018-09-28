# Copyright © Garima Kapila

import bs4, lxml.html as html, urllib2


"""
An organism is defined by its genus (taxonomic category), epithet (specific or 
second name). Additionally, a strain/variation may be specified.
ex. Genus: Saccharomyces, Epithet: Cerevisiae, Strain: S288C
	Organism('Saccharomyces', 'Cerevisiae', 'S288C')
"""
class Organism:

	def __init__(self, genus, epithet, strain=''):

		# format genus and epithet (ex. 'saccharomyces' -> 'Saccharomyces')
		self.genus = genus[0].upper() + genus[1:].lower()
		self.epithet = epithet[0].upper() + epithet[1:].lower()
		# format strain (ex. 's288c' -> 'S288C')
		self.strain = strain.upper()
		# tuple representing (genus, epithet, strain)
		self.params = (self.genus, self.epithet, self.strain)
		# formatted string of parameters
		self.info = self.genus + ' ' + self.epithet
		self.acronym = self.genus[0] + self.epithet[0]

		# update info if strain is given
		if self.strain != '':
			self.info += ', Strain ' + self.strain

		# formats organism name, ex. Saccharomyces-Cerevisiae-S288C
		self.name = '-'.join(self.params)
		if self.name[-1] == '-':
			self.name = self.name[:-1]

	"""
	returns True if there are existing proteins in the UniProt database that 
	correspond to the organism's genus, epithet, strain; False otherwise
	"""
	def is_valid(self):

		# url used for searching organism in UniProt database
		url = 'http://www.uniprot.org/uniprot/?query=' +                      \
			  '{}+{}+strain%3A{}'.format(*self.params)
		
		# read page contents of url
		soup = str(bs4.BeautifulSoup(urllib2.urlopen(url), 'html.parser'))
		tree = html.fromstring(soup)

		# return whether or not an entryID element for a protein exists
		entry_IDs = tree.xpath('//td[@class="entryID"]')
		exists = len(entry_IDs) > 0
		
		return exists

	"""
	todo, for now just hard-coded a few
	"""
	def get_taxon_ID(self):
		IDs = {
			'HS': 9606,
			'SC': 559292, # Strain S288C
			'SP': 284812, # Strain 972
			'EC': 83333, # strain K12
			'CE': 6239,
			'MM': 10090,
			'DM': 7227,
			'AT': 3702
		}
		return IDs[self.acronym]


