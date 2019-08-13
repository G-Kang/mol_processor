#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute the transition densities for excitations from the ground state
# and between excited states
# from Reimers-MOPAC INDO/SCI output
# Rebecca Gieseking, 12/16/2015
# Kevin Kang, 7/5/2018 --- minor bugs resolved
#-------------------------------------------------------------------------

import string,sys,math,re
from mol_assign import *
from constants import *
import getopt

# Initialize defaults
nstate_in = range(0,16)
ist    = 0	# Initial state
fst_in = 'auto'	# Automatically find states
wp     = [3.10] # in eV, 400nm
#wp     = [x/100.0 for x in range(200,401)]
#wp_tpa = [x/100.0 for x in range(100,401)]
erange = 1.5	# States to be included: 0 ~ 1.5X of Final-Initial energy
Te     = 200	# Entanglement time range
Ae     = 4.0 # Entanglement area, in cm^2
tstep  = 1	# Entanglement time step
#kap    = [0.0,0.1]	# Linewidth of excited states, in eV
kap = [0.0,0.1]
line   = 'l'
pol    = [-1,4,-1]

prog = ' '

helpfile = """
etpa_inter.py -i <inputfile> -o <outputfile> -n <number of states> -w <pump energy> -k <linewidth of states> -f <final state> -l <lineshape function> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -n    Number of states (including ground state)		    Default = All states
    -w    Energy of the pump laser				    Defaul  = 3.0996 eV (400 nm)
    -k    Linewidth of states (all excited states)		    Default = 0.0 eV
    -f	  Final state						    Default = 'auto'
    -l    Lineshape function (lorentzian or gaussian)		    Default = 'lorentzian' or 'l'

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:n:w:k:f:l:z:",['--help','--input=','--output=','--nstate=','--wp=','--kappa=','--final=','--line=','--program='])
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
    elif opt in ('-w','--wp'):
	wp = float(arg)
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

print 'Check the following ETPA parameters in the script'
print 'input nstate =',nstate_in,'\nwp     =',wp,'(eV)\nerange =',erange,'\nTe     = 0 ~',Te,'(fs)\ntstep  =',tstep,'(fs)\nkappa  =',kap,'(eV)\n'

[mol, program] = open_mol(infilename, prog)

if program == 'mopac':
    mol.read_atoms()
    mol.read_orbs()        # Read in the orbital coefficients and associated atoms from output
    mol.read_configs()     # Read in the configuration information
#    mol.write_dipmat(outfilename)
mol.print_states(min(nstate_in,11))
[nstate,fst] = mol.calc_nst(nstate_in,max(wp),fst_in,erange)
mol.calc_interst_tdip(nstate,program,infilename)
#mol.write_interst_tdip(nstate,infilename)

for k in kap:
    for w in wp:
	for n in nstate_in:
	    mol.calc_etpa_cs(n,w,ist,fst_in,erange,k,Te,Ae,tstep,line,pol,0.0)
	    mol.write_etpa_cs(outfilename+'_n'+str(n)+'_f'+str(fst_in[0])+'_wp'+'{0:.2f}'.format(w)+'_kap'+'{0:.2f}'.format(k))
print 'Exiting'
mol.file.close()

#-------------------------------------------------------------------------
# Using the parameters obtained above, propagate the density matrix
# with Liouville dynamics
# Liouville equation without dissipation: d rho / d t = -i/hbar [H, rho]
# v1.0 (without coupling between excited states)
#
# Kevin Kang, 7/5/2018
#-------------------------------------------------------------------------
'''
import numpy as np
import constants as cnst

# Constants
hbar = planck/(2.0*math.pi)

# Test parameters
dt = 1.0         # time step (in fs)
time = 1000.0    # total time duration for density propagation (in fs)
nstep = int(time/dt)  # number of time steps

dmat = np.zeros((nstate, nstate, nstep))
hamil = np.zeros((nstate,nstate))

dmat[0,0,0] = 1 # Density matrix at t = 0 (Only ground state)

for i in range(0,nstate)
'''
