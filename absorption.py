#!/usr/bin/python

#-------------------------------------------------------------------------
# Print stick/lorentzian abosorption spectra
# from Reimers-MOPAC INDO/SCI output
# Kevin Kang, 10/2/2018
#-------------------------------------------------------------------------

import string,sys,math,re
from mol_assign import *
from constants import *
import getopt

nstate = 'all'
emax = 10.0
emin = 0.0
gamma = 0.2

prog = ' '

helpfile = """
absorption.py -i <inputfile> -o <outputfile> -n <number of states> -z <program>

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
	nstate = int(arg)
    elif opt in ('-z','--program'):
        prog = arg.lower()

[mol, program] = open_mol(infilename, prog)

mol.exc_lrtz(nstate,emax,emin,gamma,program,infilename)
print 'Exiting'
mol.file.close()
