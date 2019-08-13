import sys
import string
import math
from properties import *
from molecule import *
from constants import *

# Holder for all properties of a system
class MopacMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        
        self.read_to('MOPAC   7.10')
        self.read_to('CHARGE','**************************************')
        self.charge = int(self.line.split()[-1]) if 'CHARGE' in self.line else 0

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('CARTESIAN COORDINATES')
        self.read_for(4)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 1:
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

    # Read the number of molecular orbitals and HOMO level
    def read_norbs(self):
        if self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('SHELL')
	self.read_for(1)
#	self.line = self.file.readline()
	line = self.line.split()
	self.homo = 0
	for i in range(1,len(line)):
	    self.homo += int(line[i])

        self.read_to('TOTAL')
        line = self.line.split()
        self.numorb = 0
        for i in range(1,len(line)):
            self.numorb += int(line[i])

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"
	print str(self.numorb),"orbitals read ( HOMO =",self.homo,")\n"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        if self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('Dipole moment=')
	line = self.line.split()
#	self.gdip = [float(line[2]),float(line[3]),float(line[4]),float(line[6])]

        # Save that dipole has been read
        if 'Dipole moment=' in self.line:
            self.dipole_read = True
            #print "Dipole read"
        else:
            print "Error: Ground-state dipole not read"

    # Read indo orbital parameters from parameter file and add hybridization terms
    # 0 = atomic number, 1 = principal quantum number n, 3 = sp zeta, 4-5 = d zetas, 6-7 = d coeffs
    def read_indoparams(self, indo = 'INDO1S.par'):
	par = open(indo,'r')
	self.orb_params = []
	for i in range(0,len(self.at_types)):
	    self.orb_params.append([])
	pline = par.readline()
	for i in range(0,80):
	    line = pline.split()
	    if line[1] in self.at_types:
		# Read data for this element
		k = self.at_types.index(line[1])
		self.orb_params[k].append(int(line[0])) # Element index
		self.orb_params[k].append(int(line[3]))	# Period number (n)
		pline = par.readline()
		line = pline.split()
		for j in range(0,5):
		    self.orb_params[k].append(float(line[j]))
		for j in range(0,6):
		    pline = par.readline()
	    else:
		for j in range(0,7):
		    pline = par.readline()
	par.close()

	self.hybrid = []
	for i in range(0,self.mol_size):
	    if self.atoms[i].elem == 'H':
        	# No hybridization
	        self.hybrid.append([])
	    else:
        	# Add sp hybridization terms
        	# Retrieve atomic orbital info
        	at_index = self.at_types.index(self.atoms[i].elem)
        	n      = self.orb_params[at_index][1]
        	spzeta = self.orb_params[at_index][2]
        	dzeta1 = self.orb_params[at_index][3]
        	dzeta2 = self.orb_params[at_index][4]
        	dcoef1 = self.orb_params[at_index][5]
        	dcoef2 = self.orb_params[at_index][6]

        	j2 = 2*n + 1
        	sp = - au2ang * j2 / (spzeta * 2.0 * math.sqrt(3.0))

        	self.hybrid.append([[sp,0.,0.],[0.,sp,0.],[0.,0.,sp]])

        	if dzeta1 > 0.0 or dzeta2 > 0.0:
		    # Add pd hybridization terms
            	    # Ordering - outer loop (px,py,pz), inner loop (z2,xy2,xy,xz,yz)
            	    j3= j2 - 2
            	    k3= j2
            
            	    t1 = math.sqrt( (2.0*spzeta)**j2 * (2.0*dzeta1)**j3 / (math.factorial(j2)*math.factorial(j3)*5.0) ) * (math.factorial(k3) / (spzeta+dzeta1)**k3 )
            	    t2 = math.sqrt( (2.0*spzeta)**j2 * (2.0*dzeta2)**j3 / (math.factorial(j2)*math.factorial(j3)*5.0) ) * (math.factorial(k3) / (spzeta+dzeta2)**k3 )
            	    pd = (t1*dcoef1 + t2*dcoef2) * au2ang

            	    # px
            	    self.hybrid[i].append([pd/math.sqrt(3.0), 0., 0.])
            	    self.hybrid[i].append([-pd, 0., 0.])
           	    self.hybrid[i].append([ 0.,-pd, 0.])
            	    self.hybrid[i].append([ 0., 0.,-pd])
            	    self.hybrid[i].append([ 0., 0., 0.])
            	    # py
             	    self.hybrid[i].append([0., pd/math.sqrt(3.0), 0.])
             	    self.hybrid[i].append([ 0., pd, 0.])
             	    self.hybrid[i].append([-pd, 0., 0.])
             	    self.hybrid[i].append([ 0., 0., 0.])
             	    self.hybrid[i].append([ 0., 0.,-pd])
             	    # pz
             	    self.hybrid[i].append([ 0., 0., -2.0*pd/math.sqrt(3.0)])
             	    self.hybrid[i].append([ 0., 0., 0.])
             	    self.hybrid[i].append([ 0., 0., 0.])
             	    self.hybrid[i].append([-pd, 0., 0.])
             	    self.hybrid[i].append([ 0.,-pd, 0.])

        if self.orb_params and self.hybrid:
            self.indoparams_read = True
        else:
            print "Error: INDO parameters not read"

    # Read the molecular orbitals
    def read_orbs(self, types=[], mo_min=1, mo_max=0):
        # Make sure atoms and norbs have been read
        if not self.atoms_read:
            self.read_atoms()
        if not self.norbs_read:
            self.read_norbs()
        if self.configs_read or self.states_read:
            self.file.seek(0)

        if mo_max == 0:
            mo_max = self.numorb

        self.mos = []
	self.orb_atom = []
        self.read_to('eigvals(Ev)')
        orb_count = 0
        aos_read = False
        orb_start = False
        orb_done = False

        # Enter main loop to get molecular orbitals
        while orb_done == False and len(self.line) > 0:
            line = self.line.split()
            orbline = len(line) - 1
            orb_count += orbline
          #  print line, orb_count, orbline

            if mo_min <= orb_count:
                # Read the relevant orbital info
		
                # Set up indices for loops later on
                if orb_start == False:
                    orb_index_a = mo_min - orb_count - 1
                else:
                    orb_index_a = -orbline
                if mo_max <= orb_count:
                    orb_index_b = mo_max - orb_count
                else:
                    orb_index_b = 0
       
                #print orb_index_a, orb_index_b
 
                # Add each MO to the list with an energy
                for i in range(1,len(line)):
		    if len(self.mos) < self.homo: occ = 2
		    else:			  occ = 0
                    self.mos.append(MolOrb(float(line[i]),int(occ)))

                self.read_for(4) 
                        
                # Loop over each AO in this set of MOs
                for i in range(0,self.numorb):
                    line = self.line.split()
                    if len(line) == 0:
                        self.line = self.file.readline()
                        line = self.line.split()
		    at = int(self.line[4:7]) - 1
		    aotype = self.line[9:12].strip()
		    if len(self.orb_atom) < self.numorb:
    		        self.orb_atom.append(at)
		    elif len(self.orb_atom) == self.numorb:
			self.orb_atom.append(self.mol_size)

		    # Store the MO coefficients
                    # WARNING - ordering of i/j is reversed from previous scripts!!! ############
		    for j in range(orb_index_a,orb_index_b):
                        self.mos[j].add_coeff(line[j])
    
                    # First orbital set only - store atomic orbital types
                    if aos_read == False:
                       # at = int(self.line[4:7]) - 1
                       # aotype = self.line[9:12].strip()
                        if 'P' in aotype: aotype = 'P'
                        if 'D' in aotype or 'xy2' in aotype: aotype = 'D' 

			if len(types) == 0: # GK: Add AO types without type file (the types of AOs are not True/False in this statement)
                            if 'S' or 'P' or 'D' in aotype:
                                self.atoms[at].add_ao(aotype)
                            else:
                                self.atoms[at].add_ao('Dxy2')

			else:
                            # Categorize all AOs by types from type file (if type file exists)
                            a = []
                            for j in types:
                                if 'Atom' in j:
                                    a.append(True) if (at >= (int(j[1])-1) and at < int(j[2])) else a.append(False)
                                elif 'Attype' in j:
                                    a.append(True) if self.atoms[at].elem in j else a.append(False)
                                elif 'Orb' in j:
                                    a.append(True) if aotype in j else a.append(False)
                                elif 'Orbtype' in j:
                                    a.append(True) if aotype in j and self.atoms[at].elem in j else a.append(False)

                            self.atoms[at].add_ao(a)
                            #print at, aotype, a

                    # Store the MO coefficients
                    # WARNING - ordering of i/j is reversed from previous scripts!!! ############
                    #for j in range(orb_index_a,orb_index_b):
                    #    self.mos[j].add_coeff(line[j])
                    self.read_for(1)

                aos_read = True
                self.read_to('eigvals(Ev)','/')
                orb_start = True
                if mo_max <= orb_count:
                    orb_done = True

            else:
                # Skip to the next section
                self.read_to('eigvals(Ev)','/')

        # Compute MO character based on AO types
        for i in self.mos:
            for m in range(0,len(types)):
                i.char.append(0.0)
                ocycle = 0
                for j in self.atoms:
                    for k in j.ao:
                        if k[m]: i.char[m] += i.coeff[ocycle]**2
                        ocycle += 1
            #print i.char

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

        if not self.orbs_read:
            self.read_orbs(types)
	
	print 'Reading Configurations...'

        self.read_to('spin-adapted configurations')
        self.read_for(4)

        self.configs = []
        line = self.line.split()
        while len(self.line) > 5:
            try:
                occ = int(self.line[53:57]) - 1
                vir = int(self.line[62:68]) - 1
            except ValueError:
                occ = 0
                vir = 0

            self.configs.append(Config(occ,vir))

	    # GK: Comment this block if type file is not needed
	    '''
            char = []
            achar = []
            bchar = []
            for i in range(0,len(self.mos[0].char)):
                char.append(self.mos[vir].char[i] - self.mos[occ].char[i])

                # achar = local character in orbitals included in type
                # bchar = local character in orbitals  outside of type
                achar.append(    min(self.mos[vir].char[i],self.mos[occ].char[i]))
                bchar.append(1 - max(self.mos[vir].char[i],self.mos[occ].char[i]))

                # RLG temporary edit for ligand clusters
                #char.append(self.mos[occ].char[i])
            self.configs.append(Config(occ,vir,char))
            self.configs[-1].achar = achar
            self.configs[-1].bchar = bchar

            #if len(self.configs) < 20:
            #    print char[0], achar[0], bchar[0], abs(char[0])+achar[0]+bchar[0]
	    '''

            self.line = self.file.readline()
            line = self.line.split()

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            print "Configurations read\n"
        else:
            print "Error: Configurations not read\n"

    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
#        if not self.configs_read and len(types) > 0:
	if not self.configs_read:
            self.read_configs(types)

        self.read_to('CI trans.')
        self.read_for(3)

        self.states = [State(0.0,0.0)]

        # Read energies and oscillator strengths
	# GK: Lines related to "types" all deleted
        line = self.line.split()
        while len(line) > 2:
            num = int(line[0])
            try:
                e = float(line[1])
                o = float(line[4])
            except ValueError:
                e = float(line[2])
                o = float(line[5])
            self.states.append(State(e,o))
            if o > 1.e-5:
                pol = [float(self.line[50:60]),float(self.line[61:71]),float(self.line[72:82])]
            else:
                pol = [0.0,0.0,0.0]
            #self.states[-1].calc_tdip(pol)
        
            self.line = self.file.readline()
            line = self.line.split()

        # Read CI coefficients
	self.read_to('State')
	nstate = 0

	while 'State' in self.line:
	    self.line = self.file.readline()

	    while 'Config' in self.line:
		line  = self.line.split()
		cnum  = int(line[1])-1
		coeff = float(line[2])
		self.states[nstate].add_coeff(cnum, coeff)
		self.line = self.file.readline()

	    self.read_for(2)
	    nstate += 1
	
	self.nstates = nstate
	
        # Save that states have been read
        if len(self.states) > 0:
            self.states_read = True
            #print "Excited states read"
        else:
            print "Error: Excited states not read"

    def calc_irj(self, i, j):
	irj = [0.0,0.0,0.0]
	for at in range(0,self.mol_size):
	    nstart = self.orb_atom.index(at)
	    nlast  = self.orb_atom.index(at+1)
	    for m in range(nstart,nlast):
		for k in range(0,3):
		    irj[k] += i.coeff[m]*j.coeff[m]*self.atoms[at].coord[k]
	    if (nlast - nstart) > 1:
		for m in range(1,4):
		    pind = m-1
		    for k in range(0,3):
			irj[k] -= i.coeff[nstart]*j.coeff[nstart+m]*self.hybrid[at][pind][k]
			irj[k] -= j.coeff[nstart]*i.coeff[nstart+m]*self.hybrid[at][pind][k]
	    if (nlast - nstart) > 4:
		for m in range(1,4):
		    for n in range(4,9):
			dind = (n-3)*3+m-1
			for k in range(0,3):
			    irj[k] -= i.coeff[nstart+n]*j.coeff[nstart+m]*self.hybrid[at][dind][k]
			    irj[k] -= j.coeff[nstart+n]*i.coeff[nstart+m]*self.hybrid[at][dind][k]

	for k in range(0,3): irj[k] *= cnst.ea2au

	return irj

    # Calculate dipole matrix :  mu_ij = < i | mu | j >  (i and j include all MOs)
    def calc_dipmat(self):
	if not self.indoparams_read:
	    self.read_indoparams()
        if not self.orbs_read:
            self.read_orbs()
	
	print 'Calculating dipole matrix in MO space...'
	self.dipmat = []
	self.effnucdip = [0.0,0.0,0.0,0.0] # Effective nuclear dipole excluding core electrons
	self.gdip = [0.0,0.0,0.0,0.0]

	for i in self.atoms:
	    for k in range(0,3):
		self.effnucdip[k] += i.valelec*i.coord[k]*cnst.ea2au

	for i in range(0,len(self.mos)):
	    self.dipmat.append([])
	    for j in range(0,len(self.mos)):
		if j >= i: self.dipmat[i].append(self.calc_irj(self.mos[i],self.mos[j]))
		else:      self.dipmat[i].append(self.dipmat[j][i])
		if i == j and i < self.homo:
		    for k in range(0,3):
			self.gdip[k] += self.dipmat[i][j][k]

	    if i % 50 == 49: print int(i+1),'/',len(self.mos)
	    elif i == len(self.mos)-1: print len(self.mos),'/',len(self.mos),'\n'
	for k in range(0,3):
	    self.gdip[k] = self.effnucdip[k] - 2*self.gdip[k]
	    self.gdip[3] += self.gdip[k]**2
	self.gdip[3] = self.gdip[3]**0.5
	
	if not len(self.dipmat) == 0:
	    self.dipmat_read = True

#	for k in range(0,3): self.gdip[k] += self.effnucdip[k]
#	print self.effnucdip
#	print self.gdip
	
