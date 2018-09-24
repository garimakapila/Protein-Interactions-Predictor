from Collector import *
from Results import *
import bs4, lxml.html as html, re, urllib2


def main():
	
	# Fetch organism transfer options
	interface_options, epithets = get_INSIDER_options()
	organisms = get_HINT_options(epithets)
	organisms = {organism.acronym: organism for organism in organisms}
	
	# Do interologs annotation transfer
	inputs = parse_arguments(interface_options, organisms)
	if inputs != None:
		for input_args in inputs:
			organism1, organism2, option, update = input_args
			collect_data(organism1, organism2, option, update)
			get_results(organism1, organism2, option)
			print '\n\n\n'



##################
# Helper Methods #
##################


# return organism1, organism2, and interface option
def parse_arguments(interface_options, organisms):
	# get list of organisms to transfer annotations
	update = False
	args = sys.argv[1:]
	if sys.argv[-1].upper() == 'UPDATE':
		update = True
		args = args[:-1]
	organisms_list, option = args_to_organisms(args, organisms)
	success = (len(organisms_list) > 0)
	if not success:
		organisms_list, option = get_valid_args(interface_options, organisms)
	if option != 'quit':
		# gather data and get results for each transfer
		inputs = []
		for organisms in organisms_list:
			org1, org2 = organisms
			if org1 != org2:
				inputs.append([org1, org2, option.upper(), update])
		return inputs

# collect data online, label interologs + interface information
# update means sequence alignment steps can be skipped
def collect_data(organism1, organism2, option, update):
	print 'Transferring protein interaction annotations, %s interfaces:\n> '  \
		'%s to %s' % (option.upper(), organism1.info, organism2.info)
	collector = Collector(organism1, organism2, option)
	return collector.run(update = update)

# cluster/classify
def get_results(organism1, organism2, option):
	print 'Analyzing Features'
	results = Results(organism1, organism2, option)
	results.run()



### Helper Methods for parse_arguments ###


# convert arguments to valid organisms
def get_valid_args(interface_options, organisms):
	message = 'Enter two species to transfer interactions from, or one'       \
	' species to transfer\nall other organism protein interactions from, or ' \
	'"arguments.txt" as input\narguments and an option for interfaces, '      \
	'examples:'                                                               \
	'\n> HS SC ALL\n> HS HQ\n> arguments.txt ALL\n\nOrganism options: \n'
	interface_options = '\n\t'.join(interface_options.keys())
	message += '\n'.join(['\t' + organism.acronym + '\t' + organism.info      \
		for organism in organisms.values()]) + '\nInterface options: \n' +    \
		'\t' + interface_options
	message += '\n\nEnter "quit" to exit\nNext time, you can directly run,' + \
		' for example, "python run.py HS SC HQ"\n'
	args = raw_input(message + '\n> ').split(' ')
	organisms_list, option = args_to_organisms(args, organisms)
	success = (len(organisms_list) > 0)
	
	while not success:
		if len(args) == 1 and args[0] == 'quit':
				return [], 'quit'
		args = raw_input('Invalid input\n> ').split(' ')
		organisms_list, option = args_to_organisms(args, organisms)
		success = (len(organisms_list) > 0)

	return organisms_list, option


# map arguments to organism options and option for interface type (ex. All)
def args_to_organisms(args, org_dict):
	organisms = []
	option = ''

	if len(args) == 2:
		# ex. python run.py arguments.txt ALL
		if args[0] == 'arguments.txt':
			if args[1].upper() == 'ALL' or args[1].upper() == 'HQ':
				option = args[1].upper()
				with open('arguments.txt', 'r') as f:
					for line in f:
						org1, org2 = line.strip('\n').upper().split(' ')
						print org1, org2
						organisms.append([org_dict[org1], org_dict[org2]])
		
		# ex. python run.py HS ALL
		elif args[0].upper() in org_dict:
			if args[1].upper() == 'ALL' or args[1].upper() == 'HQ':
				option = args[1].upper()
				org2 = org_dict[args[0].upper()]
				for org1 in org_dict.values():
					organisms.append([org1, org2])
	
	# ex. python run.py HS SC ALL
	elif len(args) == 3:
		if args[0].upper() in org_dict and args[1].upper() in org_dict and    \
			args[2].upper() == 'ALL' or args[2].upper() == 'HQ':
			option = args[2].upper()
			org1 = args[0].upper()
			org2 = args[1].upper()
			organisms.append([org_dict[org1], org_dict[org2]])
	
	return organisms, option



### Get transfer options ###


# get organisms and interface options from INSIDER database
def get_INSIDER_options():
	url = 'http://interactomeinsider.yulab.org/downloads.html'
	soup = str(bs4.BeautifulSoup(urllib2.urlopen(url), 'html.parser'))
	tree = html.fromstring(soup)
	# gather organism names from download links
	hrefs = tree.xpath('//td//a/@href')
	interface_options = {}
	epithets = {}
	for href in hrefs:
		file_name = href.split('/')[-1].replace('.txt', '').split('_')
		epithet = file_name[1]
		interface_option = ''.join([c for c in file_name[-1] if c.isupper()])
		epithets[epithet] = 1
		interface_options[interface_option] = 1
	return interface_options, epithets


# get organisms from HINT database
def get_HINT_options(epithets):
	url = 'http://hint.yulab.org/download/'
	soup = str(bs4.BeautifulSoup(urllib2.urlopen(url), 'html.parser'))
	tree = html.fromstring(soup)
	# gather organism names from download links
	hrefs = tree.xpath('//div[@class="am-u-sm-6"]//a//@href')
	hrefs = [href for href in hrefs if 'binary/hq' in href]
	# filter organisms that have predicted interface residues from INSIDER
	hrefs = [href for href in hrefs                                           \
		if any(word.lower() in href.lower() for word in epithets.keys())]
	# convert to Organism objects
	hrefs = [href.replace('/download/', '') for href in hrefs]
	hrefs = [href.replace('binary/hq', '') for href in hrefs]
	names = [href.replace('/', '') for href in hrefs]
	names = [tuple(re.findall('[A-Z][^A-Z]*', name)) for name in names]
	genus = [name[0] for name in names]
	names = [[name[1], ''.join(name[2:])] for name in names]
	epithet = [re.sub(r'[0-9]+', '', name[0]) for name in names]
	strain = [''.join([re.sub(r'[A-Z][a-z]+', '', name[0]), name[1]])
		for name in names]
	organisms = [Organism(g, e, s) for g, e, s in zip(genus, epithet, strain)]
	return organisms


if __name__ == "__main__":
	main()