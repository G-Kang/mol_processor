#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute the transition densities for excitations from the ground state
# Also compute the averge hole and electron energies and the % double
# character for each state
# from Reimers-MOPAC INDO/SCI output
# Rebecca Gieseking, 12/16/2015
# 
#-------------------------------------------------------------------------

import string,sys,math,re
from mol_assign import *
import getopt

prog = ' '

helpfile = """
out2in.py -i <inputfile> -o <outputfile> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:z:",['--help','--input=','--output=','--program='])
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
    elif opt in ('-z','--program'):
        prog = arg.lower()

out, prog = open_mol(infilename, prog)

out.write_vibs()

out.file.close()

