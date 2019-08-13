#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute INDO/CI sum-over-states polarizability
#
# Version 2.0
#
# Rebecca Gieseking, 1/6/2017
# 
#-------------------------------------------------------------------------

import string, sys
from mol_assign import *
import getopt

# Initialize defaults
omega = 0.0
gamma = 0.1088j
states = 0
prog = ' '
gamma2 = 0.0
tfilename = 'types'

helpfile = """
indo_alpha.py -i <inputfile> -o <outputfile> -w <omega> -g <gamma> -b <gamma2> -t <types> -s <# states> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -w    Energy at which polarizabilities are computed (eV)        Default = 0.0
    -g    Lifetime (broadening) of excited states (eV)              Default = 0.1088i (= 0.004 a.u.)
    -b    Lifetime of excited states of specific types (eV)         Default = 0.0 (= not used)
    -t    Types file, for used with -b. Only first type used.       Default = types
          Format:
             Atom    20               From atoms 1-20 to 21-max
             Attype  Ag               From Ag atoms to all other atoms
             Orbtype Ag D             From Ag D to all others
             Orb     D                From all D orbitals to all others

    -s    Number of states to include in SOS expression             Default = All states

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:w:g:b:t:s:z:",['--help','--input=','--output=','--omega=','--gamma=','--gamma2=','--types=','--states=','--program='])
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
        #infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-w','--omega'):
        omega = float(arg)
    elif opt in ('-g','--gamma'):
        gamma = complex(0.0,float(arg))
    elif opt in ('-b','--gamma2'):
        gamma2 = complex(0.0,float(arg))
    elif opt in ('-t','--types'):
        tfilename = arg
    elif opt in ('-s','--states'):
        states = int(arg)
    elif opt in ('-z','--program'):
        prog = arg.lower()

types = read_types(tfilename)

out, prog = open_mol(infilename, prog) 

out.read_states(types)
out.write_alpha(outfilename, omega, gamma, states, gamma2, types)

out.file.close()



