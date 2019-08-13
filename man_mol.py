import sys
import string
import math
import constants as cnst
from properties import *
from molecule import *

# Holder for all properties of a system
class ManMolecule(Molecule):

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        # Read past header
        while '!' in self.line:
            self.line = self.file.readline()

        self.read_to('ATOMS')
        self.read_for(1)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[1]),float(line[2]),float(line[3])],line[0]))
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
        self.numorb = 0
        print 'Not implemented for manual test system'

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        print 'Not implemented for manual test system'

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
        print 'Not implemented for manual test system'

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
        if not self.orbs_read:
            self.read_orbs()

        print 'Not implemented for manual test system'
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
        #if not self.configs_read:
        #    self.read_configs()
        # Read past header
        while '!' in self.line:
            self.line = self.file.readline()

        self.read_to('EXCITED STATES')
        self.read_for(2)

        self.states = [State(0.0,0.0)]

        # Read energies and oscillator strengths
        line = self.line.split()
        o = 0.0
        while len(line) > 2:
            e = float(line[0])
            mu = [float(line[1]),float(line[2]),float(line[3])]
            self.states.append(State(e,o))
            self.states[-1].calc_osc(mu)
            self.line = self.file.readline()
            line = self.line.split()
            #print self.states[-1].energy, self.states[-1].osc

        # Save that states have been read
        if len(self.states) > 1:
            self.states_read = True
            #print "Excited states read:", len(self.states) - 1
        else:
            print "Error: Excited states not read"

