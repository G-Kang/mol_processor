import sys
import string
import math
from properties import *
from molecule import *
from constants import *

# Holder for all properties of a system
class NwchMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        
	self.read_to('Charge           :')
        self.charge = int(self.line.split()[2])

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('XYZ format geometry')
        self.read_for(4)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[1]),float(line[2]),float(line[3])],line[0]))
            if line[1] not in self.at_types:
                self.at_types.append(line[0])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
	    self.mol_size = len(self.atoms)
            self.atoms_read = True
            print self.mol_size,'atoms read'
        else:
            print "Error: Atoms not read"

    # Read the number of molecular orbitals and HOMO level
    def read_norbs(self):
        self.numorb = 0

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
	if self.orbs_read or self.configs_read or self.states_read:
	    self.file.seek(0)

	self.gdip = [0.0,0.0,0.0,0.0]
	# Ground-state dipole is read to be 0 ONLY for SLR-TDDFT calculation 
	# because this is implemented in calculating static dipoles.
	'''
	self.read_to('Multipole analysis of the density')
	self.read_for(7)
	for k in range(0,3):
	    line = self.line.split()
	    self.gdip[k] = float(line[4])
	    self.read_for(1)
	for k in range(0,3):
	    self.gdip[3] += self.gdip[k]**2
	self.gdip[3] = self.gdip[3]*0.5

        # Save that dipole has been read
        if len(self.gdip) > 0:
            self.dipole_read = True
            #print "Dipole read"
        else:
            print "Error: Ground-state dipole not read"
	'''
    # Read the molecular orbitals
    def read_orbs(self, types=[], mo_min=1, mo_max=0):
        # Make sure atoms and norbs have been read
        if not self.atoms_read:
            self.read_atoms()
        if not self.norbs_read:
            self.read_norbs()
        if self.configs_read or self.states_read:
            self.file.seek(0)

        self.mos = []
        orb_count = 0

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print "Orbitals read"
        else:
            print "Error: Orbitals not read"

    # Read electron configurations for CI
    def read_configs(self,types=[]):
        if self.states_read:
            self.file.seek(0)

        self.configs = []

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"

    # Read excited-state info
    def read_states(self,types=[]):

	self.states = [State(0.0,0.0)]

	forfile = self.file.readlines()
	revfile = reversed(forfile)

	nline = len(forfile)
	nstate = 1
	for revline in revfile:
	    nline -= 1
	    if 'Convergence criterion met' in revline:
		nline += 5
		outline = forfile[nline]
		while 'Root' in outline:
		    line = outline.split()
		    num = int(line[1])
		    e = float(line[6])
		    nline += 5
		    line = forfile[nline].split()
		    o = float(line[3])
		    self.states.append(State(e,o))
		    nline += 2
		    outline = forfile[nline]
		    while 'Occ.' in outline:
			line = outline.split()
			cnum = [int(line[1])-1,int(line[5])-1]
			coeff = float(line[7])
			self.states[nstate].add_coeff(cnum, coeff)
			nline += 1
			outline = forfile[nline]
		    nline += 1
		    outline = forfile[nline]
		    nstate += 1


	'''
	self.read_to('Root')

	nstate = 1
	while 'Root' in self.line:
	    line = self.line.split()
	    num = int(line[1])
	    e = float(line[6])
	    if (nstate == 1) or (nstate == 2):
		e = e + 0.5
	    self.read_for(5)
	    o = float(self.line.split()[3])
	    self.states.append(State(e,o))
	    self.read_for(2)
	    while 'Occ.' in self.line:
		line  = self.line.split()
		cnum  = [int(line[1])-1,int(line[5])-1]
		coeff = float(line[7])
		self.states[nstate].add_coeff(cnum, coeff)
		self.line = self.file.readline()
	    self.line = self.file.readline()
	    nstate += 1
	'''
	self.nstates = len(self.states)

        if len(self.states) > 0:
            self.states_read = True
            #print "Excited states read"
        else:
            print "Error: Excited states not read"
