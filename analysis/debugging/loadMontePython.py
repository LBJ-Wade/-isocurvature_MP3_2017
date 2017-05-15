
# we only use the astropy table library
from astropy.io import ascii

def load(chainfile, paramfile):

	# read in the data with 
	data = ascii.read(chainfile, delimiter="\s")
	params = ascii.read(paramfile, delimiter="\t", format="no_header")['col1']

	# set up column names (read in from param file)
	data['col1'].name = 'acceptance'
	data['col2'].name = 'likelihood'
	for i in range(3,len(params)+3):
	    data['col' + str(i)].name = params[i-3]

	return data


def load_justparams(chainfile, paramfile):

	# read in the data with 
	data = ascii.read(chainfile, delimiter="\s")
	params = ascii.read(paramfile, delimiter="\t", format="no_header")['col1']

	# set up column names (read in from param file)
	data['col1'].name = 'acceptance'
	data['col2'].name = 'likelihood'
	for i in range(3,len(params)+3):
	    data['col' + str(i)].name = params[i-3]

	data.remove_column('acceptance')
	data.remove_column('likelihood')
	return data

### some debuggin code:
# data = load("../planckdata/r1.txt", '../planckdata/param')
# print data