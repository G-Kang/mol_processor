import sys
import string
import math
import constants as cnst
from properties   import *
from molecule     import *

# Holder for all properties of a system
class XYZMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.charge = 0.0
        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

	self.mol_size = int(self.line)
        self.read_for(2)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 0:
            line = self.line.split()
            self.atoms.append(Atom([float(line[1]),float(line[2]),float(line[3])],line[0]))
            if line[0] not in self.at_types:
                self.at_types.append(line[0])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
	    if self.mol_size != len(self.atoms):
		print "Please check the number of atoms in xyz file"
	    else:
                self.atoms_read = True
            #print "Atoms read"
        else:
            print "Error: Atoms not read x"

    # Read the number of molecular orbitals
    def read_norbs(self):
        self.numorb = 0
        print "Nothing to read for xyz format"

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        print "Nothing to read for xyz format"

        # Save that dipole has been read
        if 'Dipole moment=' in self.line:
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
        orb_count = 0
        print "Nothing to read for xyz format"

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print "Orbitals read"
        else:
            print "Error: Orbitals not read"


    # Read electron configurations for CI
    def read_configs(self):
        if self.states_read:
            self.file.seek(0)

        self.configs = []
        print "Nothing to read for xyz format"

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"


    # Read excited-state info
    def read_states(self,types=[]):

        self.states = []
	print "Nothing to read for xyz format"

        # Save that states have been read
        if len(self.states) > 1:
            self.states_read = True
            print len(self.states), "excited states read"
        else:
            print "Error: Excited states not read"
