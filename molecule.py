import sys
import string
import math
import struct
import constants as cnst
from properties import *

# Error function
def erf(x):
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

def dotprod(v1, v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

# dE in eV unit, laser freq and width in nm unit, returns a lineshape function in the unit of [s]
# Experimentally determined values: d(lambda) = 12 nm at lambda = 800 nm
def lineshape_gaus(dE=0.0, w=800.0, fwhm=12.0):
    dw = 2*fwhm/w**2/(cnst.nm2ev*cnst.ev2hz)
    dE = cnst.nm2ev*cnst.ev2hz*dE
    delta = math.sqrt(2/math.pi)*math.exp(-2*dE**2*dw**2)*dw
    return delta

def lineshape_lrtz(dE=0.0, w=800.0, fwhm=12.0):
    dw = fwhm/w**2/(cnst.nm2ev*cnst.ev2hz)
    dE = cnst.nm2ev*cnst.ev2hz*dE
    delta = (1/math.pi)*(0.5/dw)/(dE**2+(0.5/dw)**2)
    return delta

################################################################
# Parent molecule class                                        #
#                                                              #
# Use program-specific child classes (MopacMolecule, etc.) to  #
# input molecular information                                  #
#                                                              #
# Functions here handle post-processing of data and writing of #
# output files                                                 #
################################################################
class Molecule(object):
    def __init__(self, file):
        self.file = file
        self.line = self.file.readline()

        # Identify what has been read
        self.charge_read   = False
        self.atoms_read    = False
        self.norbs_read    = False
        self.dipole_read   = False
	self.indoparams_read = False
        self.orbs_read     = False
	self.ao_ovlp_read  = False
	self.ao_fock_read  = False
	self.mo_ovlp_read  = False
	self.mo_fock_read  = False
        self.configs_read  = False
        self.states_read   = False
	self.dipmat_read   = False
	self.interst_tdip_calc = False
	self.nst_calc	   = False

    ################################################
    #                                              #
    #  General read functions                      #
    #                                              #
    ################################################

    # Read to a certain message in the file
    # Use a gibberish message2 to make it useless unless specified
    def read_to(self, message, message2='qwertyuiopasdf'):
        while message not in self.line and message2 not in self.line and self.line:
            self.line = self.file.readline()

    # Read a specific number of lines in the file
    def read_for(self,lines):
        for i in range(0,lines):
            self.line = self.file.readline()

    ################################################
    #                                              #
    #  Writing simple output files                 #
    #                                              #
    ################################################

    # Center a molecule relative to the origin
    def recenter(self):
        if not self.atoms_read:
            self.read_atoms(types)

    # GK, 8/27/2018
	origin = [0.0,0.0,0.0]
	molmass = 0.0
	for a in self.atoms:
	    molmass += a.mass
	    for i in range(0,3):
		origin[i] += a.coord[i]*a.mass

	for i in range(0,3):
	    origin[i] /= molmass

	for a in self.atoms:
	    for i in range(0,3):
		a.coord[i] -= origin[i]

    # Using a template, write an input file using coordinates read in
    def write_input(self, outfilename, template):
        if not self.charge_read:
            self.read_charge()
        if not self.atoms_read:
            self.read_atoms()

        # A few default options for common outputs
        if 'def' in template:
            if 'xyz' in template:
                out = open(outfilename+'.xyz','w')
                out.write(str(len(self.atoms)) + '\n\n')
                for a in self.atoms:
                    out.write(string.rjust(a.elem,4))
                    for j in range(0,3):
                        out.write(string.rjust('%.5f'%a.coord[j],12))
                    out.write('\n')
                out.write('\n')
                out.close()
            elif 'mopac' in template:
                out = open(outfilename+'.dat','w')
                out.write('INDO XYZ RCI CHARGE='+str(self.charge)+' MAXCI=8000 WRTCI=2000 \
				\n%s INDO CIS\nINDO\n'%outfilename)
                for a in self.atoms:
                    out.write(string.rjust(a.elem,4))
                    for j in range(0,3):
                        out.write(string.rjust('%.5f'%a.coord[j],12) + ' 0')
                    out.write('\n')
                out.write('\n')
                out.close()
        # Generic output using templates
        # Template may include NAT, CHRG, TYPE, ELEM, XYZ, XXX, YYY, ZZZ as needed
        # Append suffix of template onto outfilename
        else:
            temp = open(template,'r')
            out = open(outfilename + template[template.rfind('.'):],'w')

            tline = temp.readline()
            while len(tline) > 0:
                # Single replacements
                if 'NAT' in tline:
                    tline = tline[:tline.find('NAT')] + str(len(self.atoms)) + tline[tline.find('NAT')+3:]
                if 'CHRG' in tline:
                    tline = tline[:tline.find('CHRG')] + str(self.charge) + tline[tline.find('CHRG')+4:]
                # Use line as a template for printing all atoms
                if 'ELEM' in tline or 'XXX' in tline:
                    for a in self.atoms:
                        line = tline
                        if 'ELEM' in tline:
                            line = line[:line.find('ELEM')] + string.rjust(a.elem,4) + line[line.find('ELEM')+4:]
                        if 'XYZ' in tline:
                            line = line[:line.find('XYZ')] + string.rjust('%.5f'%a.coord[0],12) + string.rjust('%.5f'%a.coord[1],12) + string.rjust('%.5f'%a.coord[2],12) + line[line.find('XYZ')+3:]
                        if 'XXX' in tline:
                            line = line[:line.find('XXX')] + string.rjust('%.5f'%a.coord[0],12) + line[line.find('XXX')+3:]
                        if 'YYY' in tline:
                            line = line[:line.find('YYY')] + string.rjust('%.5f'%a.coord[1],12) + line[line.find('YYY')+3:]
                        if 'ZZZ' in tline:
                            line = line[:line.find('ZZZ')] + string.rjust('%.5f'%a.coord[2],12) + line[line.find('ZZZ')+3:]
                        out.write(line)
                # Loop over atom types (use for ADF basis sets)
                elif 'TYPE' in tline:
                    for t in self.types:
                        line = tline
                        while 'TYPE' in tline:
                            line = line[:line.find('TYPE')] + string.rjust(t,4) + line[line.find('TYPE')+4:]
                        out.write(line)
                # Print line with edits as needed
                else:
                    out.write(tline)
                tline = temp.readline()
            out.close()

    # Write vibrational modes and other info needed for Raman calculations
    def write_vibs(self, outfilename, minfreq=500, maxfreq=2000):
        return

    def write_orbs(self, outfilename, types=[]):
        if not self.orbs_read:
            self.read_orbs(types)

        orb = open(outfilename+'.orb','w')
        orb.write("  Num    Energy")

        for j in types:
            orb.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        orb.write('\n')

        for i in range(0,len(self.mos)):
           orb.write(string.rjust(str(i+1), 5) + string.rjust("%.4f" % self.mos[i].energy, 10))
           for j in range(0,len(types)):
               orb.write(string.rjust("%.4f" % self.mos[i].char[j], 12))
           orb.write('\n')
        orb.close()

        #obin = open(outfilename+'.orb_bin','w')
        # Finish here??

    def write_omo(self, outfilename):
        if not self.orbs_read:
            self.read_orbs([])
        omo = open(outfilename+'.omo','w')
        omo.write(string.rjust(str(self.numorb),6) + '\n')

        nstart = 0
        perline = 10
	size = len(self.mos[0].coeff)
        
        while nstart < self.numorb:
            nend = nstart + perline
            if nend > self.numorb:
                nend = self.numorb
            for i in range(0,size):
                for j in range(nstart,nend):
                    omo.write(string.rjust('%.6f'%(self.mos[j].coeff[i]),10))
                omo.write('\n')
            nstart += perline

        omo.close()


    def write_osc(self, outfilename, types=[]):
        if not self.states_read:
            self.read_states(types)

        osc = open(outfilename+'.osc','w')
        osc.write("  Num    Energy     Osc Str")

        for j in types:
            osc.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        osc.write('\n')

        nstate = len(self.states)

        for j in range(1,nstate):
            osc.write(string.rjust(str(j+1),5) + string.rjust("%.4f" % self.states[j].energy,10) + string.rjust("%.5f" % self.states[j].osc,12))

            for i in range(0,len(types)):
                osc.write(string.rjust("%.4f" % self.states[j].char[i], 12))

            osc.write('\n')

        osc.close()


    # Compute the absorption spectrum
    def write_sigma(self, outfilename, types=[], max_e=8.0, e_step=0.02, gamma=0.2):
        # Make sure states have been read
        if not self.states_read:
            self.read_states(types)

        step_count = int(max_e/e_step) + 1
        sigma = []
        for i in range(0,step_count):
            sigma.append([i*e_step, 0.0])
            for j in range(0,len(types)*2):
                sigma[i].append(0.0)

        for j in self.states:
            # Add the cross-section to the appropriate bins
            if j.osc > 1e-6:
                for i in range(0,step_count):
                    current_e = sigma[i][0]
        
                    # Use Lorentzian line shape
                    phi = gamma/(((current_e - j.energy)**2 + gamma**2)*math.pi)
        
                    sigma[i][1] += phi * j.osc

                    for k in range(0,len(types)):
                        if j.char[k] > 0.0: 
                            sigma[i][2 + k*2] += phi * j.osc * j.char[k]
                        else:
                            sigma[i][3 + k*2] -= phi * j.osc * j.char[k]

        sig = open(outfilename+'.sigma','w')

        sig.write('  Energy        Tot Abs')
        for j in types:
            sig.write(string.rjust('To ' + j[0][:3] + ' ' + j[1], 15) + string.rjust('From ' + j[0][:3] + ' ' + j[1], 15))
        sig.write('\n')

        for i in range(0,step_count):
            sig.write(string.rjust("%.4f"%sigma[i][0],8)) 
            for j in sigma[i][1:]:
                sig.write(string.rjust("%.6f" % j ,15))
            sig.write('\n')
        sig.close()

    def compute_alpha(self, omega=0.0, gamma=0.1088j, states=0, pr=False, gamma2=0.0, types=[]):
        # Make sure states have been read
        #if not self.orbs_read:
        #    self.read_orbs([])
        if not self.states_read:
            self.read_states(types)

        alpha = [[[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]]
        #for i in range(0,len(types)*3):
        #    alpha.append([[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]])

        if states == 0 or states >= len(self.states):
            states = len(self.states) - 1

        g_state = gamma

        for k in range(1,states+1):
            s = self.states[k]
            #print k, s.char[0], s.achar[0],s.bchar[0],s.char[0]+s.achar[0]+s.bchar[0]

            # Compute state contribution to alpha
            #a_state = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]
            
            # Compute a state-specific gamma based on two input gamma values and excitation character
            if gamma2 != 0.0:
                if abs(s.char[0]) > 0.9999:
                    g_state = gamma2
                elif abs(s.char[0]) < 0.0001:
                    g_state = gamma
                else:
                    # Linear scaling - huge change in gamma for states with predominantly one character
                    #g_state = abs(s.char[0])*gamma2 + (1-abs(s.char[0]))*gamma

                    # Quadratic scaling
                    #g_state = gamma + abs(s.char[0])**2 * (gamma2 - gamma)

                    # Gaussian scaling
                    gscl = 4
                    g_state = (gamma + gamma2)/2 + erf(gscl*(abs(s.char[0])-0.5))/erf(gscl*0.5)*(gamma2 - gamma)/2

            for i in range(0,3):
                for j in range(0,3):
                    a_num = s.tdip[i]*s.tdip[j]
                    # Following Reimers code, convert to cm-1 here
                    a_denom_1 = (s.energy - g_state - omega) * cnst.ev2cm
                    a_denom_2 = (s.energy + g_state + omega) * cnst.ev2cm

                    alpha[0][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
                    '''
                    for k in range(0,len(types)):
                        # Try assigning each state's polarizability to the largest contribution
                        cmax = max(abs(s.char[k]),s.achar[k],s.bchar[k])
                        if abs(s.char[k]) == cmax:
                            alpha[k*3+1][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
                        elif s.achar[k] == cmax:
                            alpha[k*3+2][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
                        else:
                            alpha[k*3+3][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph

                        #alpha[k*3+1][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph * abs(s.char[k])
                        #alpha[k*3+2][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph * s.achar[k]
                        #alpha[k*3+3][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph * s.bchar[k]
                    '''


                    #a_state[i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
            #a_or = (a_state[0][0] + a_state[1][1] + a_state[2][2]) / 3.0
            #if abs(a_or) > 0.5 and pr == True:
            #if pr == True:
            #    print string.rjust(str(self.states.index(s)+1),5), string.rjust('%.4f'%abs(a_or),10), string.rjust('%.3f'%(s.energy-omega),10)
        if pr == True:
            for k in range(0,len(alpha)):
                print 'Total', string.rjust('%.4f'%abs((alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0),10), string.rjust('%.4f'%((alpha[k][0][0].real + alpha[k][1][1].real + alpha[k][2][2].real) / 3.0),10), string.rjust('%.4f'%((alpha[k][0][0].imag + alpha[k][1][1].imag + alpha[k][2][2].imag) / 3.0),10)

        # Orientationally average alpha
        #a_or = complex(0.0,0.0)
        #a_or +=  (alpha[0][0] + alpha[1][1] + alpha[2][2]) / 3.0
        #print a_or, abs(a_or)

        return alpha

    def write_alpha(self, outfilename, omega=0.0, gamma=0.1088j, states=0, gamma2=0.0, types=[]):
        alpha = self.compute_alpha(omega, gamma, states, True, gamma2, types)

        #a_or =  (alpha[0][0][0] + alpha[0][1][1] + alpha[0][2][2]) / 3.0

        al_out = open(outfilename + '.alpha', 'w')
        al_out.write('Omega = ' + '%.4f'%omega + ' eV   Gamma = ' + '%.4f'%gamma.imag + ' eV')
        if gamma2 != 0.0:
            al_out.write('Gamma2 = ' + '%.4f'%gamma2.imag + ' eV')
        al_out.write('\n\n')

        for k in range(0,len(alpha)):
            a_or =  (alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0

            al_out.write('Real polarizability (au)\n')
            for i in range(0,3):
                al_out.write(string.rjust('%.4f'%alpha[k][i][0].real,12) + string.rjust('%.4f'%alpha[k][i][1].real,12) + string.rjust('%.4f'%alpha[k][i][2].real,12) + '\n')
            al_out.write('\n')
            al_out.write('Imag polarizability (au)\n')
            for i in range(0,3):
                al_out.write(string.rjust('%.4f'%alpha[k][i][0].imag,12) + string.rjust('%.4f'%alpha[k][i][1].imag,12) + string.rjust('%.4f'%alpha[k][i][2].imag,12) + '\n')
            al_out.write('\n')

        
            al_out.write('Orientationally averaged polarizability (au)\n')
            al_out.write(string.rjust('%.4f'%a_or.real,12) + string.rjust('%.4f'%a_or.imag,12) + string.rjust('%.4f'%abs(a_or),12) + '\n\n\n')

        al_out.close()

    def print_states(self,nstates):
	if not self.states_read:
	    self.read_states()
	
	if nstates == 'all': nstates = min(nstates,len(self.states))
	
	print ' #  E (eV)'
	for i in range(0,nstates):
	    print '%s %s' % (string.rjust(str(i),2),string.rjust('%.3f'%self.states[i].energy,6))

    def write_dipmat(self, outfilename):
	if not self.dipmat_read:
	    self.calc_dipmat()
        dipout = open(outfilename+'.dipmat','w')

	dipout.write("GS Dipole: [%6.3f %6.3f %6.3f] %6.3f(Debye)\n\n" % (self.gdip[0],self.gdip[1],self.gdip[2],self.gdip[3]))
	dipout.write("  MO1  MO2    Dipole element(Debye)\n")

        for i in range(0,len(self.mos)):
	    for j in range(0,len(self.mos)):
                dipout.write(string.rjust(str(i+1),5) + string.rjust(str(j+1),5) + string.rjust("%.3f" % self.dipmat[i][j][0],9) + \
				string.rjust("%.3f" % self.dipmat[i][j][1],7) + string.rjust("%.3f" % self.dipmat[i][j][2],7) + '\n')
        dipout.close()

    def nw_read_interst_tdip(self, mpoles='*mpoles', output='*out', onlyground=False):
	if not self.states_read:
	    self.read_states()
	self.interst_tdip = [[self.gdip]]

	gtoe_tdip = []
	count = 0
	nline = 0
	with open(output,'r') as out:
	    forfile = out.readlines()
	    revfile = reversed(forfile)
	    nline = len(forfile)
	    for revline in revfile:
		nline -= 1
		if 'Convergence criterion met' in revline:
		    nline += 5
		    outline =  forfile[nline]
		    while 'Total times cpu:' not in outline:
			if 'Root' in outline and 'singlet' in outline:
			    nline += 2
			    outline = forfile[nline].split()
			    tdip = [0.0, 0.0, 0.0]
			    for k in range(0,3): tdip[k] = float(outline[2*k+3])
			    nline += 3
			    osc = float(forfile[nline].split()[3])
			    self.interst_tdip[0].append([tdip[0],tdip[1],tdip[2],osc])
			    if not onlyground: self.interst_tdip.append([[tdip[0],tdip[1],tdip[2],osc]])
# 			    if not onlyground: self.interst_tdip.append([[-tdip[0],-tdip[1],-tdip[2],osc]])
   			elif 'Number of moments' in outline:
			    nmom     = int(outline.split()[3])
			elif 'No.       roots in scan window' in outline:
			    self.nroots   = int(outline.split()[5])+1
	        	elif 'No.                total roots' in outline:
			    self.totroots = int(outline.split()[3])+1
			    break
			nline += 1
			outline = forfile[nline]
		    break
		
	'''
	with open(output,'r') as out:
	    outline = out.readline()
	    while outline != "":
		if 'Root' in outline and 'singlet' in outline:
		    out.readline()
		    line = out.readline().split()
		    tdip = [0.0, 0.0, 0.0]
		    for k in range(0,2): out.readline()
		    osc = float(out.readline().split()[3])
		    for i in range(0,3):
			tdip[i] = float(line[2*i+3])
		    self.interst_tdip[0].append([tdip[0],tdip[1],tdip[2],osc])
		    if not onlyground: self.interst_tdip.append([[-tdip[0],-tdip[1],-tdip[2],osc]])
    	        elif 'Number of moments' in outline:
		    nmom     = int(outline.split()[3])
		elif 'No.       roots in scan window' in outline:
		    self.nroots   = int(outline.split()[5])+1
	        elif 'No.                total roots' in outline:
		    self.totroots = int(outline.split()[3])+1
		    break
		outline = out.readline()
	'''
	if not onlyground:
	    fp = open(mpoles,'rb')
	    num = fp.read(8) # Read double-precision C raw binary numberas
	    i = 0
	    for sj in range(1,self.nroots):
	        for si in range(1,self.totroots):
		    self.interst_tdip[sj].append([0.0,0.0,0.0,0.0])
		    for k in range(0,3):
			self.interst_tdip[sj][si][k] = float(struct.unpack('d',fp.read(8))[0])
#			if sj <= si:
#			    self.interst_tdip[sj][si][k] = float(struct.unpack('d',fp.read(8))[0])
#			else:
#			    self.interst_tdip[sj][si][k] = -float(struct.unpack('d',fp.read(8))[0])
		    if si == sj:
			dstr = 0.0
			for k in range(0,3):
			    self.interst_tdip[sj][si][k] += self.gdip[k]
			    dstr += self.interst_tdip[sj][si][k]**2
		        self.interst_tdip[sj][si][3] = dstr**0.5
		    else:
		    	self.interst_tdip[sj][si][3] = self.calc_osc(self.interst_tdip[sj][si],self.states[si].energy,self.states[sj].energy)
		    for k in range(0,17): fp.read(8)

	if len(self.interst_tdip) > 0:
	    self.interst_tdip_calc = True
        else:
	    print "error: interstate transition dipoles not calculated"

    def calc_interst_tdip(self, nstate='all', prog='mopac', infile='', onlyground=False):
	if not self.dipole_read:
	    self.read_dipole()

	if prog=='nwchem':
	    mpoles=infile+'.mpoles'
	    output=infile+'.out'
	    self.nw_read_interst_tdip(mpoles, output, onlyground)

	else:
            # Make sure states have been read
	    if not self.states_read: self.read_states()
	    if not self.dipmat_read: self.calc_dipmat()
	
	    self.interst_tdip = [[[0.0,0.0,0.0]]]
	    self.nroots = len(self.states) if nstate == 'all' else nstate

	    # I: intial and J: final CI single excited states (I = 0 or J = 0: ground state)
	    # MU(IJ)_ia,jb = c(J)_jb * ( c(I)_ja * mu_ab - c(I)_ib * mu_ij )
	    # MU: interstate transition dipole, mu: element of the MO transition dipole matrix
	    nroots = 1 if onlyground else self.nroots

	    for sj in range(0,nroots):
	        self.interst_tdip.append([])
	        for si in range(0,self.nroots):
		    self.interst_tdip[sj].append([0.0,0.0,0.0])
		    if si == sj:
		        for k in range(0,3): self.interst_tdip[sj][si][k] += self.gdip[k]
		    for coeffj in self.states[sj].coeff:
		        j     = self.configs[coeffj[0]].occ
		        b     = self.configs[coeffj[0]].vir
		        cjb_j = coeffj[1]
		        for coeffi in self.states[si].coeff:
			    i     = self.configs[coeffi[0]].occ
			    a     = self.configs[coeffi[0]].vir
			    cja_i = self.states[si].search_config(self.configs,j,a)[1]
			    cib_i = self.states[si].search_config(self.configs,i,b)[1]
			    if si == 0 and sj == 0: continue
			    elif si == 0 and sj != 0:
			        for k in range(0,3):
			            self.interst_tdip[sj][si][k] += coeffj[1]*self.dipmat[j][b][k]*math.sqrt(2)
			    elif si != 0 and sj == 0:
			        for k in range(0,3):
				    self.interst_tdip[sj][si][k] += coeffi[1]*self.dipmat[i][a][k]*math.sqrt(2)
			    else:
			        for k in range(0,3):
				    if i == j and a != b:
				        self.interst_tdip[sj][si][k] += cjb_j*cja_i*self.dipmat[a][b][k]
				    elif i != j and a == b:
				        self.interst_tdip[sj][si][k] -= cjb_j*cib_i*self.dipmat[i][j][k]
				    elif i == j and a == b:
				        self.interst_tdip[sj][si][k] += cjb_j*cja_i*(self.dipmat[a][a][k]-self.dipmat[i][i][k])
		    if si == sj:
		        dstr = 0.0
		        for k in range(0,3): dstr += (self.interst_tdip[sj][si][k]*cnst.au2debye)**2
		        self.interst_tdip[sj][si].append(math.sqrt(dstr))
		    else:
	                self.interst_tdip[sj][si].append(self.calc_osc(self.interst_tdip[sj][si],self.states[si].energy,self.states[sj].energy))

	if len(self.interst_tdip) > 0:
	    self.interst_tdip_calc = True
        else:
	    print "Error: Interstate transition dipoles not calculated"

    def calc_osc(self, tdip, e1, e2):
	return (tdip[0]**2 + tdip[1]**2 + tdip[2]**2) * (cnst.fosc_fact*abs(e1-e2))

    def write_interst_tdip(self, nstate, outfilename):
	if not self.interst_tdip_calc:
	    self.calc_interst_tdip()
        tdipout = open(outfilename+'.tdip','w')
	tdipout.write("St.1  St.2  E2-E1(eV)        Dip (a.u.)          Osc str.\n")

        for si in range(0,min(nstate,len(self.interst_tdip))):
	    for sj in range(0,min(nstate,len(self.interst_tdip[si]))):
		tdipout.write("%s%s%s%s%s%s%s" % (str(si).rjust(4),str(sj).rjust(6),string.rjust('%.3f'%(self.states[sj].energy-self.states[si].energy),9), \
				string.rjust('%.3f'%self.interst_tdip[si][sj][0],10),string.rjust('%.3f'%self.interst_tdip[si][sj][1],8), \
				string.rjust('%.3f'%self.interst_tdip[si][sj][2],8),string.rjust('%.3f'%self.interst_tdip[si][sj][3],10)))
		if self.interst_tdip[si][sj][3] > 0.001: tdipout.write("   v")
		tdipout.write("\n")

#	tdipout.write("\nSt.1  St.2                   Sdip (Debye)         Total\n")
	
#	for si in range(0,min(nstate,len(self.interst_tdip))):
#	    if len(self.interst_tdip[si]) == 1: break
#       	    tdipout.write("%s%s         %s%s%s%s\n" % (str(si).rjust(4),str(si).rjust(6), \
#			string.rjust('%.3f'%self.interst_tdip[si][si][0],10),string.rjust('%.3f'%self.interst_tdip[si][si][1],8), \
#			string.rjust('%.3f'%self.interst_tdip[si][si][2],8),string.rjust('%.3f'%self.interst_tdip[si][si][3],10)))
	tdipout.close()

    def calc_nst(self,nstate_in='all',wp=3.0996,fst_in='auto',erange=1.5):
	if not self.states_read:
	    self.read_states()

	delE1 = self.states[1].energy - wp
	delE2 = self.states[2].energy - wp
	snum = 1

	nstate_in = len(self.states) if nstate_in=='all' else nstate_in

	if fst_in == 'auto':
	    fst = []
	    # wp is larger than the last state energy
	    if wp > self.states[nstate_in-1].energy:
	    	snum = nstate_in-1
	    	fst.append(snum)
	    	while snum-1 > 0 and abs(self.states[nstate_in-1].energy-self.states[snum-1].energy) < 1.0e-3:
		    fst.append(snum-1)
		    snum -= 1

	    # wp is smaller than the first state energy
	    elif delE1 >= 0.0:
	    	fst.append(snum)
	    	# add degenerate states
	    	while snum+1 < nstate_in and abs(self.states[1].energy-self.states[snum+1].energy) < 1.0e-3:
		    fst.append(snum+1)
		    snum += 1

	# wp is in between two state energies
	    else:
	    	while delE1 <= 0.0:
	            if delE2 > 0.0:
		    	if abs(delE1) < delE2:
			    fst.append(snum)
			    while snum-1 > 0 and abs(self.states[fst[0]].energy-self.states[snum-1].energy) < 1.0e-3:
			    	fst.append(snum-1)
			    	snum -= 1
		        else:
			    fst.append(snum+1)
			    while snum+2 < nstate_in and abs(self.states[fst[0]].energy-self.states[snum+2].energy) < 1.0e-3:
			    	fst.append(snum+2)
			    	snum += 1
		        break
		    snum += 1
	            delE1 = self.states[snum].energy - wp
	            delE2 = self.states[snum+1].energy - wp
	else: fst = fst_in

	fst.sort()
	emax = erange*self.states[max(fst)].energy

	if nstate_in == 'auto':
	    nstate = 0
	    while nstate < len(self.states) and self.states[nstate].energy < emax:
		nstate += 1
	else: nstate = nstate_in

	if len(fst) > 0:
	    self.nst_calc = True
	    return nstate,fst
        else:
	    print "Error: Number of states to be calculated not determined"

    def calc_tpa_cs(self, nstate_in='all', wp=3.0996, ist=0, fst_in='auto', erange=1.5, kap=0.0, line='l', pol=[2,2,2]):
	if self.nst_calc:
	    nstate = nstate_in
	    fst = fst_in
	else: [nstate,fst] = self.calc_nst(self,nstate_in,wp,fst_in,erange)

	F = pol[0]
	G = pol[1]
	H = pol[2]

	self.tpa_cs = 0.0
	w0     = wp*cnst.ev2hz/2.0
	kap_hz = kap*cnst.ev2hz

	S_tpa = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
	         [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
	         [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]

	if line == 'g' or line == 'gaussian':
	    delta = lineshape_gaus()
	else:
	    delta = lineshape_lrtz()

	for f in fst:
	    for n in range(1,nstate):
	        delE1 = abs(self.states[n].energy*cnst.ev2hz - w0)
		#delE2 = (self.states[n].energy - self.states[f].energy)*cnst.ev2hz + w0
		delE2 = delE1
	        if delE1 == 0: delE1 = float(1E-15)*cnst.ev2hz
		if delE2 == 0: delE2 = float(1E-15)*cnst.ev2hz

		for a in range(0,3):
	    	    for b in range(0,3):
	            	S_tpa[a][b] += (self.interst_tdip[f][n][a]*self.interst_tdip[n][ist][b]*cnst.au2nu**2) \
					/ (delE1-0.5j*kap_hz) \
					+ (self.interst_tdip[f][n][b]*self.interst_tdip[n][ist][a]*cnst.au2nu**2) \
					/ (delE2-0.5j*kap_hz)

	for a in range(0,3):
	    for b in range(0,3):
		self.tpa_cs += (F*S_tpa[a][a].real*S_tpa[b][b].real + G*S_tpa[a][b].real*S_tpa[a][b].real \
				+ H*S_tpa[a][b].real*S_tpa[b][a].real)*cnst.tpa_fact*w0**2*delta/30
		self.tpa_cs += (F*S_tpa[a][a].imag*S_tpa[b][b].imag + G*S_tpa[a][b].imag*S_tpa[a][b].imag \
				+ H*S_tpa[a][b].imag*S_tpa[b][a].imag)*cnst.tpa_fact*w0**2*delta/30

	if self.tpa_cs >= 0.0:
	    self.tpa_cs_calc = True
        else:
	    print "Error: TPA cross section not calculated"

    def write_tpa_cs(self, outfilename):
	if not self.tpa_cs_calc:
	    self.calc_tpa_cs()

        tpaout = open(outfilename+'.tpa','w')
	tpaout.write("# wp(eV) wp(nm) TPACS(cm^4s)\n")
	for e in wp:
	    tpaout.write("%s %s  %s\n" % (str(e).rjust(5),string.rjust('%.2f'%float(1240.0/e),5),string.rjust('{:.3e}'.format(self.tpa_cs),7)))
	tpaout.close()

    # Calculate the cross section for entangled two-photon absorption with a pump wavelength in eV (default: 3.0996eV (400nm))
    # istate(0) = ground state, fstate(1) = final state
    # erange: energy range of intermediate states to be included
    # (default: 1.5 times of final ~ ground states energy)
    # Te: entanglement time (default: 4000 fs), Timestep: 1fs
    # Ae: entanglement area (default: 1.0E-4 cm^2)
    def calc_etpa_cs(self, nstate_in='all', wp=3.0996, ist=0, fst_in='auto', erange=1.5, kap=0.0, Te=4000, Ae=1.0E-4, tstep=1, line='l', pol=[-1,4,-1], tau=0.0):
	if self.nst_calc:
	    nstate = nstate_in
	    fst = fst_in
	else: [nstate,fst] = self.calc_nst(self,nstate_in,wp,fst_in,erange)

	F = pol[0]
	G = pol[1]
	H = pol[2]

	Ae     = Ae/cnst.nu2cm**2
#	w0     = wp/2.0
#	kap_hz = kap
	w0     = wp*cnst.ev2hz/2.0
	kap_hz = kap*cnst.ev2hz

	self.etpa_cs_tot = []
	self.etpa_cs_real = []
	self.etpa_cs_imag = []
	self.tpa_cs = 0.0
	self.etpa_prob = []
	if tau == 0:
	    self.Te_list = range(tstep,Te+tstep,tstep)
	else:
	    self.Te_list = range(0,tau+1)

	if line == 'g' or line == 'gaussian': delta = lineshape_gaus()
	else:				      delta = lineshape_lrtz()

	S_tpa = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
	         [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
	         [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]

	for t in range(0,len(self.Te_list)):
	    if tau == 0:
		tau1 = self.Te_list[t]*1E-15
		tau2 = self.Te_list[t]*1E-15
		denom = 1.0/(30*Ae*self.Te_list[t]*1E-15)
	    else:
		tau1 = (Te-self.Te_list[t])*1E-15
		tau2 = (Te+self.Te_list[t])*1E-15
		denom = 1.0/(30*Ae*Te*1E-15)

	    self.etpa_cs_tot.append(0.0)
	    self.etpa_cs_real.append(0.0)
	    self.etpa_cs_imag.append(0.0)
	    self.etpa_prob.append(0.0)

	    S_etpa = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
		      [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
		      [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]

	    for f in fst:
		if True:
		    n = nstate_in
#		for n in range(0,nstate):
		    delE1 = abs(self.states[n].energy*cnst.ev2hz - w0)
#		    delE1 = abs(self.states[n].energy - w0)
		    delE2 = delE1
#		    delE2 = (self.states[n].energy-self.states[f].energy)*cnst.ev2hz + w0
#		    if delE1 == 0: delE1 = float(1E-15)
#		    if delE2 == 0: delE2 = float(1E-15)
		    if delE1 == 0: delE1 = float(1E-30)*cnst.ev2hz
		    if delE2 == 0: delE2 = float(1E-30)*cnst.ev2hz

		    for a in range(0,3):
		 	for b in range(0,3):
			    S_etpa[a][b] += (self.interst_tdip[f][n][a]*self.interst_tdip[n][ist][b]*cnst.au2nu**2) \
					    *(1-math.e**(-1j*tau1*delE1-0.5*tau1*kap_hz))/(delE1-0.5j*kap_hz) \
					    +(self.interst_tdip[f][n][b]*self.interst_tdip[n][ist][a]*cnst.au2nu**2) \
					    *(1-math.e**(-1j*tau2*delE2-0.5*tau2*kap_hz))/(delE2-0.5j*kap_hz)
			    if t == 0:
		        	S_tpa[a][b] += (self.interst_tdip[f][n][a]*self.interst_tdip[n][ist][b]*cnst.au2nu**2) \
					       / (delE1-0.5j*kap_hz) \
					       + (self.interst_tdip[f][n][b]*self.interst_tdip[n][ist][a]*cnst.au2nu**2) \
					       / (delE2-0.5j*kap_hz)
	    for a in range(0,3):
		for b in range(0,3):
		    self.etpa_cs_real[t] += (F*S_etpa[a][a].real*S_etpa[b][b].real + G*S_etpa[a][b].real*S_etpa[a][b].real \
					    + H*S_etpa[a][b].real*S_etpa[b][a].real)*cnst.etpa_fact*w0**2*delta*denom
		    self.etpa_cs_imag[t] += (F*S_etpa[a][a].imag*S_etpa[b][b].imag + G*S_etpa[a][b].imag*S_etpa[a][b].imag \
					    + H*S_etpa[a][b].imag*S_etpa[b][a].imag)*cnst.etpa_fact*w0**2*delta*denom
		    self.etpa_prob[t] += (F*S_tpa[a][a].real*S_tpa[b][b].real + G*S_tpa[a][b].real*S_tpa[a][b].real \
					 + H*S_tpa[a][b].real*S_tpa[b][a].real)*cnst.etpa_fact*w0**2*delta*denom
		    self.etpa_prob[t] += (F*S_tpa[a][a].imag*S_tpa[b][b].imag + G*S_tpa[a][b].imag*S_tpa[a][b].imag \
					 + H*S_tpa[a][b].imag*S_tpa[b][a].imag)*cnst.etpa_fact*w0**2*delta*denom

		    if t == 0:
			self.tpa_cs += (2*S_tpa[a][a].real*S_tpa[b][b].real + 2*S_tpa[a][b].real*S_tpa[a][b].real \
					+ 2*S_tpa[a][b].real*S_tpa[b][a].real)*cnst.tpa_fact*w0**2*delta/30
			self.tpa_cs += (2*S_tpa[a][a].imag*S_tpa[b][b].imag + 2*S_tpa[a][b].imag*S_tpa[a][b].imag \
					+ 2*S_tpa[a][b].imag*S_tpa[b][a].imag)*cnst.tpa_fact*w0**2*delta/30


	    self.etpa_cs_tot[t] = (self.etpa_cs_real[t]+self.etpa_cs_imag[t])

	if len(self.etpa_cs_tot) >= 0:
	    self.tpa_cs_calc = True
	    self.etpa_cs_calc = True
        else:
	    print "Error: ETPA cross section not calculated"

    def write_etpa_cs(self, outfilename):
	if not self.etpa_cs_calc:
	    self.calc_etpa_cs()

        etpaout = open(outfilename+'.etpa','w')
	
	etpaout.write("# Classical TPACS (cm^4s-1) = %s\n" % string.ljust('{:.3e}'.format(self.tpa_cs),7))
	for t in range(0,len(self.Te_list)):
		etpaout.write("%s  %s  %s  %s  %s\n" % (str(self.Te_list[t]).rjust(5),string.rjust('{:.3e}'.format(self.etpa_cs_tot[t]),7), \
			      string.rjust('{:.3e}'.format(self.etpa_cs_real[t]),7),string.rjust('{:.3e}'.format(self.etpa_cs_imag[t]),7), \
			      string.rjust('{:.3e}'.format(self.etpa_prob[t]),7)))
	etpaout.close()

    def exc_lrtz(self, state_n='all', emax=10.0, emin=0.0, gamma=0.2, prog='mopac', infile=''):
	if not self.states_read:
	    self.read_states()
	
	nstate = 0
	espacing = 0.01
	print nstate
	while self.states[nstate+1].energy < emax:
	    nstate += 1
	    if nstate == len(self.states)-1: break

	nstate = nstate if state_n=='all' else min(state_n,nstate)
	self.calc_interst_tdip(nstate+1,prog,infile,True)

	absint = [0.0]*int((emax-emin)/espacing+1)

	stick_out = open('exc_stick.txt','w')
	lrtz_out  = open('exc_lrtz.txt','w')
	stick_out.write('# E/eV  E/nm     f\n')
	lrtz_out.write('# E/eV     E/nm     abs/a.u.\n')

	for i in range(1,nstate+1):
	    stick_out.write('%s %s %s\n' % (string.rjust('%.2f'%self.states[i].energy,6),string.rjust('%.2f'%float(1239.84187/self.states[i].energy),7), \
			    string.rjust('%.3f'%self.interst_tdip[0][i][3],7)))
	    for j in range(0,len(absint)):
		dE         = j*espacing + emin - self.states[i].energy
		absint[j] += self.interst_tdip[0][i][3]*gamma/((dE**2+gamma**2)*math.pi)

	for i in range(0,len(absint)):
	    E = (i+0.5)*espacing + emin
	    lrtz_out.write('%s %s %s\n' % (string.rjust('%.2f'%E,6),string.rjust('%.2f'%float(1239.84187/E),7), \
			string.rjust('%.3f'%absint[i],7)))

	stick_out.close()
	lrtz_out.close()

    def singlet_lrtz(self, state_n='all', emax=5.0, emin=0.0, gamma=0.2, prog='mopac', infile=''):
	if not self.states_read:
	    self.read_states()

	nstate = 0
	espacing = 0.01
	while self.states[nstate+1].energy < emax:
	    nstate += 1
	    if nstate == len(self.states)-1: break

	nstate = nstate if state_n=='all' else min(state_n,nstate)
	self.calc_interst_tdip(nstate+1,prog,infile)

	absint = [0.0]*int((emax-emin)/espacing+1)

	stick_out = open('singlet_stick.txt','w')
	lrtz_out  = open('singlet_lrtz.txt','w')
	stick_out.write('# E/eV  E/nm     f\n')
	lrtz_out.write('# E/eV     E/nm     abs/a.u.\n')

	for si in range(0,nstate+1):
	    for sj in range(si+1,nstate+1):
		Eij = self.states[sj].energy - self.states[si].energy
		if Eij != 0.0:
		    stick_out.write('%s %s %s\n' % (string.rjust('%.2f'%Eij,6),string.rjust('%.2f'%float(1239.84187/Eij),7), \
			    string.rjust('%.3f'%self.interst_tdip[si][sj][3],7)))
 
		    for i in range(0,len(absint)):
		        dE         = i*espacing + emin - Eij
		        absint[i] += self.interst_tdip[si][sj][3]*gamma/((dE**2+gamma**2)*math.pi)

	for i in range(0,len(absint)):
	    E = (i+0.5)*espacing + emin
	    lrtz_out.write('%s %s %s\n' % (string.rjust('%.2f'%E,6),string.rjust('%.2f'%float(1239.84187/E),7), \
			string.rjust('%.3f'%absint[i],7)))
 

	print "Stick spectrum: singlet_stick.txt"
	print "Lorentzian spectrum: singlet_lrtz.txt"

	stick_out.close()
	lrtz_out.close()

    def get_mo_fock(self, init=0, final=1):
	if not self.orbs_read: self.read_orbs()

	init = self.homo-2
	final = self.homo+2

	self.mo_fock = []

	for i in range(init,final):
	    self.mo_fock.append([])
	    for j in range(init,final):
		elem = self.mos[i].energy if i == j else 0.0
		self.mo_fock[i-init].append(elem)
		print elem,
	    print

    def rot_matrix(self, init=0, final=1, ind1=0, ind2=1, gam=0.0):
	U_dim = final-init
	U = []

	for i in range(U_dim):
	    U.append([])
	    for j in range(U_dim):
		elem = 0.0
		if   i == ind1 and j == ind1: elem = math.cos(gam)
		elif i == ind1 and j == ind2: elem = math.sin(gam)
		elif i == ind2 and j == ind1: elem = -math.sin(gam)
		elif i == ind2 and j == ind2: elem = math.cos(gam)
		U[i].append(elem)

	return U

    def rotate_mo(self, mo1, mo2, gam):
	new_mo1 = mo1
	new_mo2 = mo2

	for i in range(len(mo1.coeff)):
	    new_mo1.coeff[i] =  math.cos(gam)*mo1.coeff[i] + math.sin(gam)*mol2.coeff[i]
	    new_mo2.coeff[i] = -math.sin(gam)*mo1.coeff[i] + math.cos(gam)*mol2.coeff[i]

	return new_mo1, new_mo2
    '''
    def calc_dG(self, mo_start, mo_end):

	for s in range(mo_start-1,mo_end):
	    for t in range(i+1,mo_end):
		Ast = 0.0
		Bst = 0.0
	    	
		for k in range(0,3):
    '''
