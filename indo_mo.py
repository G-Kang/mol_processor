#!/usr/bin/python

#-------------------------------------------------------------------------
# Extract INDO MO coefficients from Mopac to be read into another Mopac job
#
# Version 2.0
#
# Rebecca Gieseking, 12/2/2016
# 
#-------------------------------------------------------------------------

import string, sys
from molecule import *
import getopt

helpfile = """
indo_mo.py -i <inputfile> -o <outputfile> 

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:",['--help','--input=','--output='])
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
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg

# Open input and output files
try:
    out = Molecule(infilename+'.log')
except IOError:
    out = Molecule(infilename+'.out')

out.write_omo(outfilename)


out.file.close()

