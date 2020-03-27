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
from operator import itemgetter
import getopt

def calc_percentage(matrix,size,f):
    outfile_sort = open('state_portion_f'+str(f)+'.out','w')
    outfile_matrix = open('tpa_matrix_f'+str(f)+'.out','w')
    row = [j for sub in matrix for j in sub]
    total = sum(row)

    indices,row_sorted = zip(*sorted(enumerate(row),reverse=True,key=itemgetter(1)))
    percentage = [i*100.0/total for i in row_sorted]

    for n in range(len(row_sorted)):
	i = indices[n]/size
	j = indices[n]%size
	tpacs = row_sorted[n]
	outfile_sort.write('i: %s, j: %s, TPACS: %s (%s)\n' % (str(i),str(j),string.rjust('{:.3e}'.format(tpacs),7),string.rjust('%.2f'%percentage[n],5)))
    for i in range(size):
	for j in range(size):
	    outfile_matrix.write('%s\t%s\t%s\n' % (str(i),str(j),string.rjust('{:.3e}'.format(matrix[i][j]),7)))

# Initialize defaults
ist    = 0	# Initial state
nst_in = 11
fst_in = [5]	# Automatically find states
wp     = [3.10] # in eV, 400nm
#wp     = [x/100.0 for x in range(200,401)]
#wp_tpa = [x/100.0 for x in range(100,401)]
erange = 1.5	# States to be included: 0 ~ 1.5X of Final-Initial energy
Te     = 200	# Entanglement time range
Ae     = 1.0E-6 # Entanglement area, in cm^2
tstep  = 1	# Entanglement time step
kap    = [0.0,0.1]	# Linewidth of excited states, in eV
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
	wp = [float(arg)]
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
print 'input nstate =',nst_in,'\nwp     =',wp,'(eV)\nerange =',erange,'\nTe     = 0 ~',Te,'(fs)\ntstep  =',tstep,'(fs)\nkappa  =',kap,'(eV)\n'

[mol, program] = open_mol(infilename, prog)

if program == 'mopac':
    mol.read_atoms()
    mol.read_orbs()        # Read in the orbital coefficients and associated atoms from output
    mol.read_configs()     # Read in the configuration information
#    mol.write_dipmat(outfilename)
mol.print_states(11)
mol.calc_interst_tdip(16,program,infilename)
mol.write_interst_tdip(16,infilename)

tpa_matrix = [[0.0 for x in range(nst_in)] for y in range(nst_in)]
for f in fst_in:
    for k in kap:
	for w in wp:
	    for i in range(nst_in):
		for j in range(nst_in):
		    etpafile = outfilename+'_n'+str(i)+'-'+str(j)+'_f'+str(f)+'_wp'+'{0:.2f}'.format(w)+'_kap'+'{0:.2f}'.format(k)
		    etpa_matrix = mol.calc_etpa_sr([i,j],w,ist,f,erange,k,Te,Ae,tstep,line,pol,0.0)
		    tpa_matrix[i][j] = etpa_matrix[5]
		    mol.write_etpa_cs(etpa_matrix,etpafile)
    calc_percentage(tpa_matrix,nst_in,f)
print 'Exiting'
mol.file.close()


