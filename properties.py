import math
import constants as cnst

# Atom
class Atom(object):
    def __init__(self, coord, elem, at_chg=0.0):
        self.coord = coord
        self.elem = elem
	self.at_chg = at_chg
        self.ao = []
#	self.mo = []
        self.mass = cnst.mass_dict[self.elem]
	self.valelec = cnst.valelec_dict[self.elem]

    # Add atomic orbitals
    def add_ao(self, type):
        self.ao.append(type)

#    def add_mo(self,mo = MolOrb()):
#	self.mo.append(mo)

# Molecular orbital
class MolOrb(object):
    def __init__(self, energy=0.0, occ=0.0):
        self.energy = energy
        self.coeff = []
	self.orb_atom = []
        self.occ = occ
        self.char = []

    # Add more AO coefficients
    def add_coeff(self, new_coeff):
        self.coeff.append(float(new_coeff))

    def add_orb_atom(self, at):
	self.orb_atom.append(int(at))

# Electron configuration
class Config(object):
    def __init__(self, occ, vir, char = []):
#    def __init__(sefl, occ, vir, char):
        self.change_occ(occ)
        self.change_vir(vir)
        self.char = char

    # Change occupied orbital(s) this configuration excites from
    def change_occ(self, occ):
	self.occ = occ
	'''
        if occ is list:
            self.occ = occ
        else:
            self.occ = [occ]
	'''

    # Change virtual orbital(s) this configuration excites to
    def change_vir(self, vir):
	self.vir = vir
	'''
        if vir is list:
            self.vir = vir
        else:
            self.vir = [vir]
	'''

# Excited state
class State(object):
    def __init__(self, energy=0.0, osc=0.0):
        self.energy = energy
        self.osc = osc
        self.coeff = []
        self.char = []
        self.tdip = [0.0,0.0,0.0]

    def add_coeff(self, config, coeff):
        self.coeff.append([config, coeff])

    def calc_tdip(self, pol):
        for i in range(0,3):
            self.tdip[i] = pol[i]*math.sqrt(self.osc/(cnst.fosc_fact*self.energy))
    '''
    def calc_osc(self, tdip):
        self.tdip = tdip
        self.osc = (tdip[0]**2 + tdip[1]**2 + tdip[2]**2) * (cnst.fosc_fact*self.energy)
        #print self.energy, self.tdip
    '''

    def search_config(self, configs, occ, vir):
        for i in self.coeff:
	    if configs[i[0]].occ == occ and configs[i[0]].vir == vir:
		return i
	return [0, 0.0]

