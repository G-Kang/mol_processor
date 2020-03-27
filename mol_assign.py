from mopac_mol import *
from adf_mol   import *
from qchem_mol import *
from xyz_mol   import *
from nwch_mol  import *
from gaus_mol  import *
from man_mol   import *

def open_mol(filename, prog = ' '):
    # Open input and output files
    try:
        file = open(filename,'r')
    except IOError:
        try:
            file = open(filename+'.out','r')
        except IOError:
	    try:
                file = open(filename+'.log','r')
	    except:
		file = open(filename+'.xyz','r')
 
    # Check what program this input is associated with and open as the appropriate class
    if prog == ' ':
	if file.readline().rstrip('\n').isdigit() and len(file.readline().rstrip('\n')) == 1:
	    prog = 'xyz'
	    file.seek(0)
	else:
            fline = file.read(2000).lower()
            file.seek(0)
            if 'mopac: public domain version' in fline:
                prog = 'mopac'
            elif 'amsterdam density functional' in fline or 'amsterdam modeling suite' in fline:
                prog = 'adf'
            elif 'q-chem, inc., pleasanton, ca' in fline:
                prog = 'qchem'
	    elif 'northwest computational chemistry' in fline:
		prog = 'nwchem'
	    elif 'gaussian' in fline:
		prog = 'gaussian'
            else:
                # Assume manual input
                prog = 'man'

    if prog == 'mopac':
        out = MopacMolecule(file)
    elif prog == 'adf':
        out = AdfMolecule(file)
    elif prog == 'qchem':
        out = QchemMolecule(file)
    elif prog == 'nwchem':
	out = NwchMolecule(file)
    elif prog == 'gaussian':
	out = GausMolecule(file)
    elif prog == 'xyz':
	out = XYZMolecule(file)
    elif prog == 'man':
        out = ManMolecule(file)
    else:
        sys.exit("Program "+prog+" not recognized")

    return out, prog

def read_types(tfilename):
    types = []
    try:
        tfile = open(tfilename,'r')
        for line in tfile.readlines():
            print 'Type   %s' % line,
            types.append(line.split())
    except IOError:
        print '%s not found; no types assigned' % tfilename

    return types


