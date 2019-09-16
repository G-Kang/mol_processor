#!/usr/bin/python

import string,sys,math,re
from mol_assign import *
from constants import *
import getopt

# Initialize defaults
nstate_in = 'all'# Number of excited states
ist    = 0	# Initial state
#fst_in = 'auto'	# Final state
fst_in = [1,2,6]
wp     = [x/100.0 for x in range(100,401)] # in eV
erange = 1.5	# States to be included: 0 ~ 1.5X of Final-Initial energy
kap    = [0.0,0.1]
line   = 'l'
pol    = [2,2,2]

prog = ' '

helpfile = """
tpa.py -i <inputfile> -o <outputfile> -n <number of states> -k <linewidth of states> -f <final state> -l <lineshape func> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -n    Number of states (including ground state)		    Default = All states
    -k	  Linewidth of states (all excited states)		    Default = 0.0 eV
    -f	  Final state						    Default = 'auto'
    -l	  Lineshape function (lorentzian or gaussian)		    Default = 'lorentzian' or 'l'

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:n:k:f:l:z:",['--help','--input=','--output=','--nstate=','--kappa=','--final=','--line=','--program='])
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
    elif opt in ('-k','--kappa'):
	kap = float(arg)
    elif opt in ('-f','--final'):
	try: fst_in = [int(arg)]
	except ValueError:
	    if ',' in arg:
	        fst_in = arg.split(',')
	    elif arg == 'auto':
	        fst_in = str(arg)
    elif opt in ('-l','--line'):
	line = str(arg)
    elif opt in ('-z','--program'):
        prog = arg.lower()

print 'Check the following TPA parameters in the script'
print 'input nstate =',nstate_in,'\nwp     =',wp[0],'~',wp[-1],'(eV)\nerange =',erange,'\nkappa =',kap,'(eV)\n'

[mol, program] = open_mol(infilename, prog)

if program == 'mopac':
    mol.read_atoms()
    mol.read_orbs()        # Read in the orbital coefficients and associated atoms from output
    mol.read_configs()     # Read in the configuration information
mol.print_states(min(nstate_in,11))
[nstate,fst] = mol.calc_nst(nstate_in,max(wp),fst_in,erange)
mol.calc_interst_tdip(nstate,program,infilename)
mol.write_interst_tdip(nstate,infilename)

for f in fst_in:
    for k in kap:
	tpaout = open(outfilename+'_f'+str(f)+'_kap'+'{0:.2f}'.format(k)+'.tpa','w')
	tpaout.write("# wp(eV) wp(nm) TPACS(cm^4s)\n")
	for w in wp:
	    [nstate,fst] = mol.calc_nst(nstate,w,[f],erange)
	    mol.calc_tpa_cs(nstate,w,ist,fst,erange,k,line,pol)
	    tpaout.write("%s %s  %s\n" % (str(w).rjust(5),string.rjust('%.2f'%float(1240.0/w),5),string.rjust('{:.3e}'.format(mol.tpa_cs),7)))
	tpaout.close()

print 'Exiting'
mol.file.close()
