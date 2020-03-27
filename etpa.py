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
nstate_in = 'all'
ist    = 0	# Initial state
#fst_in = 'auto'	# Automatically find states
#fst_in = [i+1 for i in range(10)]
fst_in = range(8,9)
#fst_in = [2,3,5]
#wp     = [3.10] # in eV, 400nm
wp     = [x/100.0 for x in range(150,351)]
#wp	= [x/100.0 for x in range(500,901)]
#wp_tpa	= [x/100.0 for x in range(500,901)]

erange = 10.0	# States to be included: 0 ~ 1.5X of Final-Initial energy
Te     = 200	# Entanglement time range
Ae     = 1.0E-6 # Entanglement area, in cm^2
tstep  = 1	# Entanglement time step
#kap    = [0.0,0.1]	# Linewidth of excited states, in eV
kap = [0.0,0.1]
line   = 'l'
pol    = [-1,4,-1]

prog = ' '

helpfile = """
etpa.py -i <inputfile> -o <outputfile> -n <number of states> -w <pump energy> -k <linewidth of states> -f <final state> -l <lineshape function> -z <program>

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
mol.write_interst_tdip(nstate,infilename)

for f in fst_in:
    for k in kap:
	tpaout = open(outfilename+'_f'+str(f)+'_kap'+'{0:.2f}'.format(k)+'.tpa','w')
	tpaout.write("# wp(eV) wp(nm) TPACS(cm^s)\n")

	etpaout = open(outfilename+'_f'+str(f)+'_kap'+'{0:.2f}'.format(k)+'.etpa','w')
	spmout = open(outfilename+'_f'+str(f)+'_kap'+'{0:.2f}'.format(k)+'.spm','w')

	etpa_matrix = []
	spm_matrix = []
	for w in wp:
	    [nstate,fst] = mol.calc_nst(nstate_in,w,[f],erange)
	    etpaval = mol.calc_etpa_cs(nstate,w,ist,fst,erange,k,Te,Ae,tstep,line,pol,0.0)
	    etpa_matrix.append(etpaval[1])
	    spm_matrix.append(etpaval[4])
#	    mol.write_etpa_cs(outfilename+'_f'+str(f)+'_wp'+'{0:.2f}'.format(w)+'_kap'+'{0:.2f}'.format(k))
	    tpaout.write("%s %s  %s\n" % (str(w).rjust(5),string.rjust('%2f'%float(1240.00/w),5),string.rjust('{:.3e}'.format(etpaval[5]),7)))
	tpaout.close()
	for t in etpaval[0]:
	    for w in wp:
		i = etpaval[0].index(t)
		j = wp.index(w)
		etpaout.write("  %s" % str(etpa_matrix[j][i]))
		spmout.write("  %s" % str(spm_matrix[j][i]))
	    etpaout.write("\n")
	    spmout.write("\n")
	etpaout.close()
	spmout.close()
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
