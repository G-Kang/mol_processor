#!/usr/bin/python

#-------------------------------------------------------------------------
# Calculate spontaneous emission lifetime
# Kevin Kang, 9/27/2018
#-------------------------------------------------------------------------

import string,sys,math,re
from mol_assign import *
from constants import *
import getopt

fstate = 1

prog = ' '

helpfile = """
em_tau.py -i <inputfile> -o <outputfile> -f <final state> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -f    Final excited state					    Default = 1st excited state

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:f:z:",['--help','--input=','--output=','--fstate=','--program='])
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
    elif opt in ('-f','--fstate'):
	fstate = int(arg)
    elif opt in ('-z','--program'):
        prog = arg.lower()

[mol, program] = open_mol(infilename, prog)

for f in range(1,fstate+1):
        print 'f=',f,' ',mol.calc_em_tau(0,f,program,infilename),'Hz'
mol.file.close()
