#!/usr/bin/python

import string,sys,math,re
from mol_assign import *
from constants import *
import getopt

# Initialize defaults
nstate_in = 'all'# Number of excited states
prog = ' '

helpfile = """
tdip.py -i <inputfile> -o <outputfile> -n <number of states> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -n    Number of states (including ground state)		    Default = All states

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:n:z:",['--help','--input=','--output=','--nstate=','--program='])
except getopt.GetoptError as err:
    print str(err)
    print helpfile
    sys.exit(2)

if remainder != []:
    print "Error: options not read - %r" % remainder
    print helpfile
    sys.exit(2)


for opt, arg in options:
    if opt in ('-h','--help'):
        print helpfile
        sys.exit(2)
    elif opt in ('-i','--input'):
        infilename = arg
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-n','--nstate'):
	nstate_in = int(arg)
    elif opt in ('-z','--program'):
        prog = arg.lower()

print 'input nstate =',nstate_in,'\n'

[mol, program] = open_mol(infilename, prog)

if program == 'mopac':
    mol.read_atoms()
    mol.read_orbs()        # Read in the orbital coefficients and associated atoms from output
    mol.read_configs()     # Read in the configuration information
mol.print_states(min(nstate_in,11))
mol.calc_interst_tdip(nstate_in,program,infilename)
mol.write_interst_tdip(nstate_in,infilename)

print 'Exiting'
mol.file.close()
