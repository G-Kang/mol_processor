import sys
import string
import math
import constants as cnst
from properties import *
from molecule import *

# Holder for all properties of a system
class QchemMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('$molecule')
        self.line = self.file.readline()
        self.charge = int(self.line.split()[0])

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"


    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        # Check whether the output is an optimization (read optimized geometry) or not (read initial geometry)
        self.read_to('$molecule')
	self.read_to('$rem')
        while 'end' not in self.line and len(self.line) > 0:
            opt = True if 'opt' in self.line.lower() else False
            self.line = self.file.readline()

        self.atoms = []
        self.at_types = []

        if opt:
            self.read_to("OPTIMIZATION CONVERGED")
        self.read_to('Coordinates (Angstroms)','Standard Nuclear Orientation')
        self.read_for(2) if 'Coordinates (Angstroms)' in self.line else self.read_for(3)

        while '-----------' not in self.line and len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[2]),float(line[3]),float(line[4])],line[1]))
            if line[1] not in self.at_types:
                self.at_types.append(line[1])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
	    self.mol_size = len(self.atoms)
            self.atoms_read = True
            print self.mol_size,'atoms read'
        else:
            print "Error: Atoms not read"

    # Read the number of molecular orbitals
    def read_norbs(self):
        if self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('Nuclear Repulsion Energy')
	self.read_for(1)
	self.homo = int(self.line.split()[2])

        self.read_for(2)
        line = self.line.split()
        self.numao = int(self.line.split()[5])

        # Save that norbs has been read
        if self.numao > 0:
            self.numao_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"
	print str(self.numao),"AOs read ( HOMO =",self.homo,")"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        if self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

	self.read_to('Dipole Moment (Debye)')
	self.read_for(2)
	self.gdip = self.line.split()[1]

        # Save that dipole has been read
        if self.gdip:
            self.dipole_read = True
            #print "Dipole read"
        else:
            print "Error: Ground-state dipole not read"

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
	self.orb_atom = []
	orb_count = 0
	aos_read = False
	orb_done = False
	self.read_to('MOLECULAR ORBITAL COEFFICIENTS')

	while not orb_done:
	    self.read_to('eigenvalues:')
	    line = self.line.split()
	    orbline = len(line) - 1
	    orb_count += orbline

	    for i in range(1,len(line)):
		occ = 2 if len(self.mos) < self.homo else 0
		self.mos.append(MolOrb(float(line[i]),int(occ)))
	    self.read_for(1)
	    line = self.line.split()
	
	    if not aos_read:
	        at = 0
	        oldline = line
	        elem = line[1]

	    for i in range(0,self.numao):
		if not aos_read:
		    if not isinstance(line[2], int): line.insert(2,1)
		    if not (line[1] == oldline[1] and line[2] == oldline[2]): at += 1
		    self.orb_atom.append(at)

		    aotype = line[3]
		    if   'p' in aotype: aotype = 'p'
		    elif 'd' in aotype: aotype = 'd'
		    self.atoms[at].add_ao(aotype)
		    oldline = line

		elif len(self.orb_atom) == self.numao:
		    self.orb_atom.append(self.mol_size)

		for j in range(-orbline,0):
		    self.mos[j].add_coeff(line[j])

		line = self.file.readline().split()
	    aos_read = True

	    if len(line) == 0: orb_done = True

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
	    self.numorb = len(self.mos)
            print self.numorb,"MOs read"
        else:
            print "Error: Orbitals not read"

    def read_ao_matrix(self, overlap=True, fock=True):
        if not self.atoms_read:
            self.read_atoms()
        if not self.orbs_read:
            self.read_orbs()

	self.file.seek(0)

	if overlap:
	    self.read_to('Overlap Matrix')
	    self.ao_ovlp = []
	    for i in range(0,self.numao): self.ao_ovlp.append([])
	    ao_count = 0
	    while ao_count < self.numao:
	        self.read_for(1)
	        line = self.line.split()
	        aoline = len(line)
	        ao_count += aoline
		    
	        for i in range(0,self.numao):
		    self.read_for(1)
		    line = self.line.split()
		    for j in range(-aoline,0):
		        self.ao_ovlp[i].append(float(line[j]))

            if self.ao_ovlp > 0:
                self.ao_ovlp_read = True
                print "AO Overlap Matrix read"
            else:
                print "Error: AO Overlap Matrix not read"

	if fock:
	    self.read_to('Fock Matrix')
	    self.ao_fock = []
	    for i in range(0,self.numao): self.ao_fock.append([])
	    ao_count = 0
	    while ao_count < self.numao:
	        self.read_for(1)
	        line = self.line.split()
	        aoline = len(line)
	        ao_count += aoline
		    
	        for i in range(0,self.numao):
		    self.read_for(1)
		    line = self.line.split()
		    for j in range(-aoline,0):
		        self.ao_fock[i].append(float(line[j]))

            if self.ao_fock > 0:
                self.ao_fock_read = True
                print "AO Fock Matrix read"
            else:
                print "Error: AO Fock Matrix not read"

    def calc_mo_matrix(self, init=0, final=1, overlap=True, fock=True):
	if (overlap and not self.ao_ovlp_read) or (fock and not self.ao_fock_read):
	    self.read_ao_matrix(overlap,fock)

	init = self.homo-2
	final = self.homo+2
	
	if overlap: self.mo_ovlp = []
	if fock:    self.mo_fock = []
	
	for i,imo in enumerate(self.mos[init:final]):
	    if overlap: self.mo_ovlp.append([])
	    if fock:    self.mo_fock.append([])
	    for j,jmo in enumerate(self.mos[init:final]):
		if overlap: self.mo_ovlp[i].append(0.0)
		if fock:    self.mo_fock[i].append(0.0)
		if j >= i:
		    for m in range(0,len(imo.coeff)):
			if overlap: self.mo_ovlp[i][j] += imo.coeff[m]**2*self.ao_ovlp[m][m]
			if fock:    self.mo_fock[i][j] += imo.coeff[m]**2*self.ao_fock[m][m]
		        for n in range(0,m):
			    if overlap: self.mo_ovlp[i][j] += 2*imo.coeff[m]*jmo.coeff[n]*self.ao_ovlp[m][n]
			    if fock:    self.mo_fock[i][j] += 2*imo.coeff[m]*jmo.coeff[n]*self.ao_fock[m][n]
		else:
		    if overlap: self.mo_ovlp[i][j] = self.mo_ovlp[j][i]
		    if fock:    self.mo_fock[i][j] = self.mo_fock[j][i]
		print i,j,self.mo_ovlp[i][j]

    # Read electron configurations for CI
    def read_configs(self):
        if self.states_read:
            self.file.seek(0)

        if not self.orbs_read:
            self.read_orbs()

        print 'Not yet implemented for Qchem output'
        self.configs = []

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"


    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
        if not self.configs_read:
            self.read_configs()

        print 'Not yet implemented for Qchem output'
        self.states = [State(0.0,0.0)]

        # Save that states have been read
        if len(self.states) > 0:
            self.states_read = True
            #print "Excited states read"
        else:
            print "Error: Excited states not read"

    def write_vibs(self):
        if not self.atoms_read:
            self.read_atoms()

        self.modes = []
        nmodes = 0
        
        nmo  = open(     'nmodes.inp', 'w')
        while 'STANDARD THERMODYNAMIC' not in self.line and len(self.line) > 0:
            self.read_to('Frequency:')
            nmo.write('            ' + self.line[12:] + '           --------------------   --------------------   --------------------\n')
            self.read_for(7)
            while 'TransDip' not in self.line and len(self.line) > 0:
                nmo.write(self.line)
                self.line = self.file.readline()
            self.read_for(2)
            nmo.write('\n\n')
        nmo.close()

        atm  = open( 'atomicmass.inp', 'w')
        for i in range(0,len(self.atoms)):
            atm.write(string.rjust(str(i),5) + '. ' + string.ljust(self.atoms[i].elem,3) + string.rjust('%.8f'%self.atoms[i].mass,12) + '\n')
        atm.close()




