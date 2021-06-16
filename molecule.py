import sys, string, math, struct, constants as cnst
from properties import *

def erf(x):
    sign = 1 if x >= 0 else -1
    x = abs(x)
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * math.exp(-x * x)
    return sign * y


def dotprod(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def lineshape_gaus(dE=0.0, w=800.0, fwhm=12.0):
    dw = 2 * fwhm / w ** 2 / (cnst.nm2ev * cnst.ev2hz)
    dE = cnst.nm2ev * cnst.ev2hz * dE
    delta = math.sqrt(2 / math.pi) * math.exp(-2 * dE ** 2 * dw ** 2) * dw
    return delta


def lineshape_lrtz(E=cnst.nmev / 400.0, w=800.0, fwhm=12.0):
    w_ev = cnst.nmev / w
    fwhm_ev = cnst.nmev * (1 / (w - fwhm / 2) - 1 / (w + fwhm / 2))
    dE = abs(2 * w_ev - E)
    delta = 1 / math.pi * fwhm_ev
    delta = delta / (dE ** 2 + fwhm_ev ** 2) / cnst.ev2hz
    return delta


class Molecule(object):

    def __init__(self, file):
        self.file = file
        self.line = self.file.readline()
        self.charge_read = False
        self.atoms_read = False
        self.norbs_read = False
        self.dipole_read = False
        self.indoparams_read = False
        self.orbs_read = False
        self.ao_ovlp_read = False
        self.ao_fock_read = False
        self.mo_ovlp_read = False
        self.mo_fock_read = False
        self.configs_read = False
        self.states_read = False
        self.dipmat_read = False
        self.interst_tdip_calc = False
        self.nst_calc = False
        self.tpa_cs_calc = False
        self.etpa_cs_calc = False

    def read_to(self, message, message2='qwertyuiopasdf'):
        while message not in self.line and message2 not in self.line and self.line:
            self.line = self.file.readline()

    def read_for(self, lines):
        for i in range(0, lines):
            self.line = self.file.readline()

    def recenter(self):
        if not self.atoms_read:
            self.read_atoms(types)
        origin = [
         0.0, 0.0, 0.0]
        molmass = 0.0
        for a in self.atoms:
            molmass += a.mass
            for i in range(0, 3):
                origin[i] += a.coord[i] * a.mass

        for i in range(0, 3):
            origin[i] /= molmass

        for a in self.atoms:
            for i in range(0, 3):
                a.coord[i] -= origin[i]

    def write_input(self, outfilename, template):
        if not self.charge_read:
            self.read_charge()
        if not self.atoms_read:
            self.read_atoms()
        if 'def' in template:
            if 'xyz' in template:
                out = open(outfilename + '.xyz', 'w')
                out.write(str(len(self.atoms)) + '\n\n')
                for a in self.atoms:
                    out.write(string.rjust(a.elem, 4))
                    for j in range(0, 3):
                        out.write(string.rjust('%.5f' % a.coord[j], 12))

                    out.write('\n')

                out.write('\n')
                out.close()
            elif 'mopac' in template:
                out = open(outfilename + '.dat', 'w')
                out.write('INDO XYZ RCI CHARGE=' + str(self.charge) + ' MAXCI=8000 WRTCI=2000 \t\t\t\t\n%s INDO CIS\nINDO\n' % outfilename)
                for a in self.atoms:
                    out.write(string.rjust(a.elem, 4))
                    for j in range(0, 3):
                        out.write(string.rjust('%.5f' % a.coord[j], 12) + ' 0')

                    out.write('\n')

                out.write('\n')
                out.close()
        else:
            temp = open(template, 'r')
            out = open(outfilename + template[template.rfind('.'):], 'w')
            tline = temp.readline()
            while len(tline) > 0:
                if 'NAT' in tline:
                    tline = tline[:tline.find('NAT')] + str(len(self.atoms)) + tline[tline.find('NAT') + 3:]
                if 'CHRG' in tline:
                    tline = tline[:tline.find('CHRG')] + str(self.charge) + tline[tline.find('CHRG') + 4:]
                if 'ELEM' in tline or 'XXX' in tline:
                    for a in self.atoms:
                        line = tline
                        if 'ELEM' in tline:
                            line = line[:line.find('ELEM')] + string.rjust(a.elem, 4) + line[line.find('ELEM') + 4:]
                        if 'XYZ' in tline:
                            line = line[:line.find('XYZ')] + string.rjust('%.5f' % a.coord[0], 12) + string.rjust('%.5f' % a.coord[1], 12) + string.rjust('%.5f' % a.coord[2], 12) + line[line.find('XYZ') + 3:]
                        if 'XXX' in tline:
                            line = line[:line.find('XXX')] + string.rjust('%.5f' % a.coord[0], 12) + line[line.find('XXX') + 3:]
                        if 'YYY' in tline:
                            line = line[:line.find('YYY')] + string.rjust('%.5f' % a.coord[1], 12) + line[line.find('YYY') + 3:]
                        if 'ZZZ' in tline:
                            line = line[:line.find('ZZZ')] + string.rjust('%.5f' % a.coord[2], 12) + line[line.find('ZZZ') + 3:]
                        out.write(line)

                else:
                    if 'TYPE' in tline:
                        for t in self.types:
                            line = tline
                            while 'TYPE' in tline:
                                line = line[:line.find('TYPE')] + string.rjust(t, 4) + line[line.find('TYPE') + 4:]

                            out.write(line)

                    else:
                        out.write(tline)
                tline = temp.readline()

            out.close()

    def write_vibs(self, outfilename, minfreq=500, maxfreq=2000):
        pass

    def write_orbs(self, outfilename, types=[]):
        if not self.orbs_read:
            self.read_orbs(types)
        orb = open(outfilename + '.orb', 'w')
        orb.write('  Num    Energy')
        for j in types:
            orb.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        orb.write('\n')
        for i in range(0, len(self.mos)):
            orb.write(string.rjust(str(i + 1), 5) + string.rjust('%.4f' % self.mos[i].energy, 10))
            for j in range(0, len(types)):
                orb.write(string.rjust('%.4f' % self.mos[i].char[j], 12))

            orb.write('\n')

        orb.close()

    def write_omo(self, outfilename):
        if not self.orbs_read:
            self.read_orbs([])
        omo = open(outfilename + '.omo', 'w')
        omo.write(string.rjust(str(self.numorb), 6) + '\n')
        nstart = 0
        perline = 10
        size = len(self.mos[0].coeff)
        while nstart < self.numorb:
            nend = nstart + perline
            if nend > self.numorb:
                nend = self.numorb
            for i in range(0, size):
                for j in range(nstart, nend):
                    omo.write(string.rjust('%.6f' % self.mos[j].coeff[i], 10))

                omo.write('\n')

            nstart += perline

        omo.close()

    def write_osc(self, outfilename, types=[]):
        if not self.states_read:
            self.read_states(types)
        osc = open(outfilename + '.osc', 'w')
        osc.write('  Num    Energy     Osc Str')
        for j in types:
            osc.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        osc.write('\n')
        nstate = len(self.states)
        for j in range(1, nstate):
            osc.write(string.rjust(str(j + 1), 5) + string.rjust('%.4f' % self.states[j].energy, 10) + string.rjust('%.5f' % self.states[j].osc, 12))
            for i in range(0, len(types)):
                osc.write(string.rjust('%.4f' % self.states[j].char[i], 12))

            osc.write('\n')

        osc.close()

    def write_sigma(self, outfilename, types=[], max_e=8.0, e_step=0.02, gamma=0.2):
        if not self.states_read:
            self.read_states(types)
        step_count = int(max_e / e_step) + 1
        sigma = []
        for i in range(0, step_count):
            sigma.append([i * e_step, 0.0])
            for j in range(0, len(types) * 2):
                sigma[i].append(0.0)

        for j in self.states:
            if j.osc > 1e-06:
                for i in range(0, step_count):
                    current_e = sigma[i][0]
                    phi = gamma / (((current_e - j.energy) ** 2 + gamma ** 2) * math.pi)
                    sigma[i][1] += phi * j.osc
                    for k in range(0, len(types)):
                        if j.char[k] > 0.0:
                            sigma[i][(2 + k * 2)] += phi * j.osc * j.char[k]
                        else:
                            sigma[i][(3 + k * 2)] -= phi * j.osc * j.char[k]

        sig = open(outfilename + '.sigma', 'w')
        sig.write('  Energy        Tot Abs')
        for j in types:
            sig.write(string.rjust('To ' + j[0][:3] + ' ' + j[1], 15) + string.rjust('From ' + j[0][:3] + ' ' + j[1], 15))

        sig.write('\n')
        for i in range(0, step_count):
            sig.write(string.rjust('%.4f' % sigma[i][0], 8))
            for j in sigma[i][1:]:
                sig.write(string.rjust('%.6f' % j, 15))

            sig.write('\n')

        sig.close()

    def compute_alpha(self, omega=0.0, gamma=complex(0.0, 0.1088), states=0, pr=False, gamma2=0.0, types=[]):
        if not self.states_read:
            self.read_states(types)
        alpha = [
         [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
          [
           complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
          [
           complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]]
        if states == 0 or states >= len(self.states):
            states = len(self.states) - 1
        g_state = gamma
        for k in range(1, states + 1):
            s = self.states[k]
            if gamma2 != 0.0:
                if abs(s.char[0]) > 0.9999:
                    g_state = gamma2
                elif abs(s.char[0]) < 0.0001:
                    g_state = gamma
                else:
                    gscl = 4
                    g_state = (gamma + gamma2) / 2 + erf(gscl * (abs(s.char[0]) - 0.5)) / erf(gscl * 0.5) * (gamma2 - gamma) / 2
            for i in range(0, 3):
                for j in range(0, 3):
                    a_num = s.tdip[i] * s.tdip[j]
                    a_denom_1 = (s.energy - g_state - omega) * cnst.ev2cm
                    a_denom_2 = (s.energy + g_state + omega) * cnst.ev2cm
                    alpha[0][i][j] += (a_num / a_denom_1 + a_num / a_denom_2) * cnst.alph

        if pr == True:
            for k in range(0, len(alpha)):
                print 'Total', string.rjust('%.4f' % abs((alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0), 10), string.rjust('%.4f' % ((alpha[k][0][0].real + alpha[k][1][1].real + alpha[k][2][2].real) / 3.0), 10), string.rjust('%.4f' % ((alpha[k][0][0].imag + alpha[k][1][1].imag + alpha[k][2][2].imag) / 3.0), 10)

        return alpha

    def write_alpha(self, outfilename, omega=0.0, gamma=complex(0.0, 0.1088), states=0, gamma2=0.0, types=[]):
        alpha = self.compute_alpha(omega, gamma, states, True, gamma2, types)
        al_out = open(outfilename + '.alpha', 'w')
        al_out.write('Omega = ' + '%.4f' % omega + ' eV   Gamma = ' + '%.4f' % gamma.imag + ' eV')
        if gamma2 != 0.0:
            al_out.write('Gamma2 = ' + '%.4f' % gamma2.imag + ' eV')
        al_out.write('\n\n')
        for k in range(0, len(alpha)):
            a_or = (alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0
            al_out.write('Real polarizability (au)\n')
            for i in range(0, 3):
                al_out.write(string.rjust('%.4f' % alpha[k][i][0].real, 12) + string.rjust('%.4f' % alpha[k][i][1].real, 12) + string.rjust('%.4f' % alpha[k][i][2].real, 12) + '\n')

            al_out.write('\n')
            al_out.write('Imag polarizability (au)\n')
            for i in range(0, 3):
                al_out.write(string.rjust('%.4f' % alpha[k][i][0].imag, 12) + string.rjust('%.4f' % alpha[k][i][1].imag, 12) + string.rjust('%.4f' % alpha[k][i][2].imag, 12) + '\n')

            al_out.write('\n')
            al_out.write('Orientationally averaged polarizability (au)\n')
            al_out.write(string.rjust('%.4f' % a_or.real, 12) + string.rjust('%.4f' % a_or.imag, 12) + string.rjust('%.4f' % abs(a_or), 12) + '\n\n\n')

        al_out.close()

    def print_states(self, nstates):
        if not self.states_read:
            self.read_states()
        if nstates == 'all':
            nstates = min(nstates, len(self.states))
        print ' #  E (eV)'
        for i in range(0, nstates):
            print '%s %s' % (string.rjust(str(i), 2), string.rjust('%.3f' % self.states[i].energy, 6))

    def write_dipmat(self, outfilename):
        if not self.dipmat_read:
            self.calc_dipmat()
        dipout = open(outfilename + '.dipmat', 'w')
        dipout.write('GS Dipole: [%6.3f %6.3f %6.3f] %6.3f(Debye)\n\n' % (self.gdip[0], self.gdip[1], self.gdip[2], self.gdip[3]))
        dipout.write('  MO1  MO2    Dipole element(Debye)\n')
        for i in range(0, len(self.mos)):
            for j in range(0, len(self.mos)):
                dipout.write(string.rjust(str(i + 1), 5) + string.rjust(str(j + 1), 5) + string.rjust('%.3f' % self.dipmat[i][j][0], 9) + string.rjust('%.3f' % self.dipmat[i][j][1], 7) + string.rjust('%.3f' % self.dipmat[i][j][2], 7) + '\n')

        dipout.close()

    def nw_read_interst_tdip(self, mpoles='*mpoles', output='*out', onlyground=False):
        if not self.states_read:
            self.read_states()
        self.interst_tdip = [[self.gdip]]
        gtoe_tdip = []
        count = 0
        nline = 0
        with open(output, 'r') as (out):
            forfile = out.readlines()
            revfile = reversed(forfile)
            nline = len(forfile)
            for revline in revfile:
                nline -= 1
                if 'Convergence criterion met' in revline:
                    nline += 5
                    outline = forfile[nline]
                    while 'Total times' not in outline:
                        if 'Root' in outline and 'singlet' in outline:
                            nline += 2
                            outline = forfile[nline].split()
                            tdip = [0.0, 0.0, 0.0]
                            for k in range(0, 3):
                                tdip[k] = float(outline[(2 * k + 3)])

                            nline += 3
                            osc = float(forfile[nline].split()[3])
                            self.interst_tdip[0].append([tdip[0], tdip[1], tdip[2], osc])
                            if not onlyground:
                                self.interst_tdip.append([[tdip[0], tdip[1], tdip[2], osc]])
                        else:
                            if 'Number of moments' in outline:
                                nmom = int(outline.split()[3])
                            else:
                                if 'No.       roots in scan window' in outline:
                                    self.nroots = int(outline.split()[5]) + 1
                                else:
                                    if 'No.                total roots' in outline:
                                        self.totroots = int(outline.split()[3]) + 1
                                        break
                        nline += 1
                        outline = forfile[nline]

                    break

        if not onlyground:
            fp = open(mpoles, 'rb')
            num = fp.read(8)
            i = 0
            for sj in range(1, self.nroots):
                for si in range(1, self.totroots):
                    self.interst_tdip[sj].append([0.0, 0.0, 0.0, 0.0])
                    for k in range(0, 3):
                        self.interst_tdip[sj][si][k] = float(struct.unpack('d', fp.read(8))[0])

                    if si == sj:
                        dstr = 0.0
                        for k in range(0, 3):
                            self.interst_tdip[sj][si][k] += self.gdip[k]
                            dstr += self.interst_tdip[sj][si][k] ** 2

                        self.interst_tdip[sj][si][3] = dstr ** 0.5
                    else:
                        self.interst_tdip[sj][si][3] = self.calc_osc(self.interst_tdip[sj][si], self.states[si].energy, self.states[sj].energy)
                    for k in range(0, 17):
                        fp.read(8)

        if len(self.interst_tdip) > 0:
            self.interst_tdip_calc = True
        else:
            print 'error: interstate transition dipoles not calculated'

    def calc_interst_tdip(self, nstate='all', prog='mopac', infile='', onlyground=False):
        if not self.dipole_read:
            self.read_dipole()
        if prog == 'nwchem':
            mpoles = infile + '.mpoles'
            output = infile + '.out'
            self.nw_read_interst_tdip(mpoles, output, onlyground)

	if prog == 'adf':
	   if not self.states_read:
		self.read_states() 

        if prog == 'mopac':
            if not self.states_read:
                self.read_states()
            if not self.dipmat_read:
                self.calc_dipmat()
            self.interst_tdip = [[[0.0, 0.0, 0.0]]]
            self.nroots = len(self.states) if nstate == 'all' else nstate
            nroots = 1 if onlyground else self.nroots
            for sj in range(0, nroots):
                self.interst_tdip.append([])
                for si in range(0, self.nroots):
                    self.interst_tdip[sj].append([0.0, 0.0, 0.0])
                    if si == sj:
                        for k in range(0, 3):
                            self.interst_tdip[sj][si][k] += self.gdip[k]

                    for coeffj in self.states[sj].coeff:
                        j = self.configs[coeffj[0]].occ
                        b = self.configs[coeffj[0]].vir
                        cjb_j = coeffj[1]
                        for coeffi in self.states[si].coeff:
                            i = self.configs[coeffi[0]].occ
                            a = self.configs[coeffi[0]].vir
                            cja_i = self.states[si].search_config(self.configs, j, a)[1]
                            cib_i = self.states[si].search_config(self.configs, i, b)[1]
                            if si == 0 and sj == 0:
                                continue
                            elif si == 0 and sj != 0:
                                for k in range(0, 3):
                                    self.interst_tdip[sj][si][k] += coeffj[1] * self.dipmat[j][b][k] * math.sqrt(2)

                            elif si != 0 and sj == 0:
                                for k in range(0, 3):
                                    self.interst_tdip[sj][si][k] += coeffi[1] * self.dipmat[i][a][k] * math.sqrt(2)

                            else:
                                for k in range(0, 3):
                                    if i == j and a != b:
                                        self.interst_tdip[sj][si][k] += cjb_j * cja_i * self.dipmat[a][b][k]
                                    elif i != j and a == b:
                                        self.interst_tdip[sj][si][k] -= cjb_j * cib_i * self.dipmat[i][j][k]
                                    elif i == j and a == b:
                                        self.interst_tdip[sj][si][k] += cjb_j * cja_i * (self.dipmat[a][a][k] - self.dipmat[i][i][k])

                    if si == sj:
                        dstr = 0.0
                        for k in range(0, 3):
                            dstr += (self.interst_tdip[sj][si][k] * cnst.au2debye) ** 2

                        self.interst_tdip[sj][si].append(math.sqrt(dstr))
                    else:
                        self.interst_tdip[sj][si].append(self.calc_osc(self.interst_tdip[sj][si], self.states[si].energy, self.states[sj].energy))

        if len(self.interst_tdip) > 0:
            self.interst_tdip_calc = True
        else:
            print 'Error: Interstate transition dipoles not calculated'

    def calc_osc(self, tdip, e1, e2):
        return (tdip[0] ** 2 + tdip[1] ** 2 + tdip[2] ** 2) * (cnst.fosc_fact * abs(e1 - e2))

    def write_interst_tdip(self, nstate, outfilename):
        if not self.interst_tdip_calc:
            self.calc_interst_tdip()
        tdipout = open(outfilename + '.tdip', 'w')
        tdipout.write('St.1  St.2  E2-E1(eV)        Dip (a.u.)          Osc str.\n')
        for si in range(0, min(nstate, len(self.interst_tdip))):
            for sj in range(0, min(nstate, len(self.interst_tdip[si]))):
                tdipout.write('%s%s%s%s%s%s%s' % (str(si).rjust(4), str(sj).rjust(6), string.rjust('%.3f' % (self.states[sj].energy - self.states[si].energy), 9),
                 string.rjust('%.3f' % self.interst_tdip[si][sj][0], 10), string.rjust('%.3f' % self.interst_tdip[si][sj][1], 8),
                 string.rjust('%.3f' % self.interst_tdip[si][sj][2], 8), string.rjust('%.3f' % self.interst_tdip[si][sj][3], 10)))
                if self.interst_tdip[si][sj][3] > 0.001:
                    tdipout.write('   v')
                tdipout.write('\n')

        tdipout.close()

    def calc_nst(self, nstate_in='all', wp=3.0996, fst_in='auto', erange=1.5):
        if not self.states_read:
            self.read_states()
        delE1 = self.states[1].energy - wp
        delE2 = self.states[2].energy - wp
        snum = 1
        nstate_in = len(self.states) if nstate_in == 'all' else nstate_in
        if fst_in == 'auto':
            fst = []
            if wp > self.states[(nstate_in - 1)].energy:
                snum = nstate_in - 1
                fst.append(snum)
                while snum - 1 > 0 and abs(self.states[(nstate_in - 1)].energy - self.states[(snum - 1)].energy) < 0.001:
                    fst.append(snum - 1)
                    snum -= 1

            elif delE1 >= 0.0:
                fst.append(snum)
                while snum + 1 < nstate_in and abs(self.states[1].energy - self.states[(snum + 1)].energy) < 0.001:
                    fst.append(snum + 1)
                    snum += 1

            else:
                while delE1 <= 0.0:
                    if delE2 > 0.0:
                        if abs(delE1) < delE2:
                            fst.append(snum)
                            while snum - 1 > 0 and abs(self.states[fst[0]].energy - self.states[(snum - 1)].energy) < 0.001:
                                fst.append(snum - 1)
                                snum -= 1

                        else:
                            fst.append(snum + 1)
                            while snum + 2 < nstate_in and abs(self.states[fst[0]].energy - self.states[(snum + 2)].energy) < 0.001:
                                fst.append(snum + 2)
                                snum += 1

                        break
                    snum += 1
                    delE1 = self.states[snum].energy - wp
                    delE2 = self.states[(snum + 1)].energy - wp

        else:
            fst = fst_in
        fst.sort()
        emax = erange * self.states[max(fst)].energy
        if nstate_in == 'auto':
            nstate = 0
            while nstate < len(self.states) and self.states[nstate].energy < emax:
                nstate += 1

        else:
            nstate = nstate_in
        if len(fst) > 0:
            self.nst_calc = True
            return (
             nstate, fst)
        print 'Error: Number of states to be calculated not determined'

    def calc_tpa_cs(self, nstate_in='all', wp=3.0996, ist=0, fst_in='auto', erange=1.5, kap=0.0, line='l', pol=[2, 2, 2]):
        if self.nst_calc:
            nstate = nstate_in
            fst = fst_in
        else:
            nstate, fst = self.calc_nst(self, nstate_in, wp, fst_in, erange)
        F = pol[0]
        G = pol[1]
        H = pol[2]
        w0 = wp / 2.0
        self.tpa_cs = 0.0
        S_tpa = [
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
        delta = 1 / (0.1 * math.pi)
        for f in fst:
            for n in range(1, nstate):
                delE1 = abs(self.states[n].energy / cnst.au2ev - w0)
                delE2 = delE1
                if delE1 == 0:
                    delE1 = float(1e-30)
                if delE2 == 0:
                    delE2 = float(1e-30)
                if f == n:
                    tdip1 = [ self.interst_tdip[f][n][k] - self.interst_tdip[ist][ist][k] for k in range(0, 3) ]
                else:
                    tdip1 = self.interst_tdip[f][n]
                tdip2 = self.interst_tdip[n][ist]
                for a in range(0, 3):
                    for b in range(0, 3):
                        S_tpa[a][b] += tdip1[a] * tdip2[b] / (delE1 - complex(0.0, 0.5) * kap) + tdip1[b] * tdip2[a] / (delE2 - complex(0.0, 0.5) * kap)

        for a in range(0, 3):
            for b in range(0, 3):
                self.tpa_cs += (F * S_tpa[a][a].real * S_tpa[b][b].real + G * S_tpa[a][b].real * S_tpa[a][b].real + H * S_tpa[a][b].real * S_tpa[b][a].real) * cnst.tpa_fact * w0 ** 2 * delta / 30

        if self.tpa_cs >= 0.0:
            self.tpa_cs_calc = True
        else:
            print 'Error: TPA cross section not calculated'

    def write_tpa_cs(self, outfilename):
        if not self.tpa_cs_calc:
            self.calc_tpa_cs()
        tpaout = open(outfilename + '.tpa', 'w')
        tpaout.write('# wp(eV) wp(nm) TPACS(cm^4s)\n')
        for e in wp:
            tpaout.write('%s %s  %s\n' % (str(e).rjust(5), string.rjust('%.2f' % float(1240.0 / e), 5), string.rjust(('{:.3e}').format(self.tpa_cs), 7)))

        tpaout.close()

    def calc_etpa_cs(self, nstate_in='all', wp=3.0996, ist=0, fst_in='auto', erange=1.5, kap=0.0, Te=4000, Ae=1e-06, tstep=1, line='l', pol=[-1, 4, -1], tau=0.0):
        if self.nst_calc:
            nstate = nstate_in
            fst = fst_in
        else:
            nstate, fst = self.calc_nst(self, nstate_in, wp, fst_in, erange)
        F = pol[0]
        G = pol[1]
        H = pol[2]
        w0 = wp / 2.0
        etpa_cs_tot = []
        etpa_cs_real = []
        etpa_cs_imag = []
        tpa_cs = 0.0
        etpa_prob = []
        if tau == 0:
            Te_list = range(tstep, Te + tstep, tstep)
        else:
            Te_list = range(-tau, tau + 1)
        delta_tpa = 1 / (0.1 * math.pi)
        S_tpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                 [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                 [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
        for t in range(0, len(Te_list)):
            if tau == 0:
                tau1 = Te_list[t] * 1e-15
                tau2 = Te_list[t] * 1e-15
                denom = 1.0 / (30 * Ae * Te_list[t] * 1e-15)
            else:
                tau1 = (Te - Te_list[t]) * 1e-15
                tau2 = (Te + Te_list[t]) * 1e-15
                denom = 1.0 / (30 * Ae * Te * 1e-15)
            etpa_cs_tot.append(0.0)
            etpa_cs_real.append(0.0)
            etpa_cs_imag.append(0.0)
            etpa_prob.append(0.0)
            S_etpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                      [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                      [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
            for f in fst:
                delta = cnst.ev2hz / self.calc_em_tau(0,f)
                for n in range(0, nstate):
                    delE1 = abs(self.states[n].energy - w0)
                    delE2 = abs(self.states[n].energy - w0)
                    if delE1 == 0:
                    #    continue
                        delE1 = float(1e-30)
                    if delE2 == 0:
                        delE2 = float(1e-30)
                    if f == n:
                        tdip1 = [ self.interst_tdip[f][n][k] - self.interst_tdip[ist][ist][k] for k in range(0, 3) ]
                    else:
                        tdip1 = self.interst_tdip[f][n]
                    tdip2 = self.interst_tdip[n][ist]
                    for a in range(0, 3):
                        for b in range(0, 3):
                            S_etpa[a][b] += tdip1[a] * tdip2[b] * (1 - math.e ** ((complex(0.0, -1.0) * tau1 * delE1 - 0.5 * tau1 * kap) * cnst.ev2hz)) / (delE1 - complex(0.0, 0.5) * kap) + tdip1[b] * tdip2[a] * (1 - math.e ** ((complex(0.0, -1.0) * tau2 * delE2 - 0.5 * tau2 * kap) * cnst.ev2hz)) / (delE2 - complex(0.0, 0.5) * kap)
                            if t == 0:
                                S_tpa[a][b] += tdip1[a] * tdip2[b] / (delE1 - complex(0.0, 0.5) * kap) + tdip1[b] * tdip2[a] / (delE2 - complex(0.0, 0.5) * kap)

            for a in range(0, 3):
                for b in range(0, 3):
                    etpa_cs_real[t] += (F * S_etpa[a][a].real * S_etpa[b][b].real + G * S_etpa[a][b].real * S_etpa[a][b].real + H * S_etpa[a][b].real * S_etpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_cs_imag[t] += (F * S_etpa[a][a].imag * S_etpa[b][b].imag + G * S_etpa[a][b].imag * S_etpa[a][b].imag + H * S_etpa[a][b].imag * S_etpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[t] += (F * S_tpa[a][a].real * S_tpa[b][b].real + G * S_tpa[a][b].real * S_tpa[a][b].real + H * S_tpa[a][b].real * S_tpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[t] += (F * S_tpa[a][a].imag * S_tpa[b][b].imag + G * S_tpa[a][b].imag * S_tpa[a][b].imag + H * S_tpa[a][b].imag * S_tpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    if t == 0:
                        tpa_cs += (2 * S_tpa[a][a].real * S_tpa[b][b].real + 2 * S_tpa[a][b].real * S_tpa[a][b].real + 2 * S_tpa[a][b].real * S_tpa[b][a].real) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30
                        tpa_cs += (2 * S_tpa[a][a].imag * S_tpa[b][b].imag + 2 * S_tpa[a][b].imag * S_tpa[a][b].imag + 2 * S_tpa[a][b].imag * S_tpa[b][a].imag) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30

            etpa_cs_tot[t] = etpa_cs_real[t]

        if len(etpa_cs_tot) >= 0:
            self.tpa_cs_calc = True
            self.etpa_cs_calc = True
        else:
            print 'Error: ETPA cross section not calculated'
        return (
         Te_list, etpa_cs_tot, etpa_cs_real, etpa_cs_imag, etpa_prob, tpa_cs)

    def calc_etpa_sr(self, nstate_in=[
 0, 0], wp=3.0996, ist=0, fst_in=1, erange=1.5, kap=0.0, Te=4000, Ae=1e-06, tstep=1, line='l', pol=[-1, 4, -1], tau=0.0):
        i = nstate_in[0]
        j = nstate_in[1]
        f = fst_in
        F = pol[0]
        G = pol[1]
        H = pol[2]
        w0 = wp / 2.0
        etpa_cs_tot = []
        etpa_cs_real = []
        etpa_cs_imag = []
        etpa_prob = []
        tpa_cs = 0.0
        if tau == 0:
            Te_list = range(tstep, Te + tstep, tstep)
        else:
            Te_list = range(0, tau + 1)
        delta_tpa = 1 / (0.1 * math.pi)
        bi_tpa = [
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
        bj_tpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
         [
          complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]

	tau_inv = self.calc_em_tau(0,f)
	if tau_inv == 0: delta = 0.0
	else: delta = cnst.ev2hz / tau_inv

        delEi = abs(self.states[i].energy - w0)
        delEj = abs(self.states[j].energy - w0)
        if delEi == 0:
            delEi = float(1e-30)
        if delEj == 0:
            delEj = float(1e-30)
        if f == i:
            tdipi1 = [ self.interst_tdip[f][i][k] - self.interst_tdip[ist][ist][k] for k in range(0, 3) ]
        else:
            tdipi1 = self.interst_tdip[f][i]
        tdipi2 = self.interst_tdip[i][ist]
        if f == j:
            tdipj1 = [ self.interst_tdip[f][j][k] - self.interst_tdip[ist][ist][k] for k in range(0, 3) ]
        else:
            tdipj1 = self.interst_tdip[f][j]
        tdipj2 = self.interst_tdip[j][ist]
        for t in range(0, len(Te_list)):
            if tau == 0:
                tau1 = Te_list[t] * 1e-15
                tau2 = Te_list[t] * 1e-15
                denom = 1.0 / (30 * Ae * Te_list[t] * 1e-15)
            else:
                tau1 = (Te - Te_list[t]) * 1e-15
                tau2 = (Te + Te_list[t]) * 1e-15
                denom = 1.0 / (30 * Ae * Te * 1e-15)
            etpa_cs_tot.append(0.0)
            etpa_cs_real.append(0.0)
            etpa_cs_imag.append(0.0)
            etpa_prob.append(0.0)
            bi_etpa = [
             [
              complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
             [
              complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
             [
              complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
            bj_etpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
             [
              complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
             [
              complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
            for a in range(0, 3):
                for b in range(0, 3):
                    bi_etpa[a][b] = 2 * (tdipi1[a] * tdipi2[b]) * (1 - math.e ** ((complex(0.0, -1.0) * tau1 * delEi - 0.5 * tau1 * kap) * cnst.ev2hz)) / (delEi - complex(0.0, 0.5) * kap)
                    bj_etpa[a][b] = 2 * (tdipj1[a] * tdipj2[b]) * (1 - math.e ** ((complex(0.0, -1.0) * tau1 * delEj - 0.5 * tau1 * kap) * cnst.ev2hz)) / (delEj - complex(0.0, 0.5) * kap)
                    if t == 0:
                        bi_tpa[a][b] = 2 * (tdipi1[a] * tdipi2[b]) / (delEi - complex(0.0, 0.5) * kap)
                        bj_tpa[a][b] = 2 * (tdipj1[a] * tdipj2[b]) / (delEj - complex(0.0, 0.5) * kap)

            for a in range(0, 3):
                for b in range(0, 3):
                    etpa_cs_real[t] += (F * bi_etpa[a][a].real * bj_etpa[b][b].real + G * bi_etpa[a][b].real * bj_etpa[a][b].real + H * bi_etpa[a][b].real * bj_etpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_cs_imag[t] += (F * bi_etpa[a][a].imag * bj_etpa[b][b].imag + G * bi_etpa[a][b].imag * bj_etpa[a][b].imag + H * bi_etpa[a][b].imag * bj_etpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[t] += (F * bi_tpa[a][a].real * bj_tpa[b][b].real + G * bi_tpa[a][b].real * bj_tpa[a][b].real + H * bi_tpa[a][b].real * bj_tpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[t] += (F * bi_tpa[a][a].imag * bj_tpa[b][b].imag + G * bi_tpa[a][b].imag * bj_tpa[a][b].imag + H * bi_tpa[a][b].imag * bj_tpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    if t == 0:
                        tpa_cs += (2 * bi_tpa[a][a].real * bj_tpa[b][b].real + 2 * bi_tpa[a][b].real * bj_tpa[a][b].real + 2 * bi_tpa[a][b].real * bj_tpa[b][a].real) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30
                        tpa_cs += (2 * bi_tpa[a][a].imag * bj_tpa[b][b].imag + 2 * bi_tpa[a][b].imag * bj_tpa[a][b].imag + 2 * bi_tpa[a][b].imag * bj_tpa[b][a].imag) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30

            etpa_cs_tot[t] = etpa_cs_real[t]

        if len(etpa_cs_tot) >= 0:
            self.tpa_cs_calc = True
            self.etpa_cs_calc = True
        else:
            print 'Error: ETPA cross section not calculated'
        return (Te_list, etpa_cs_tot, etpa_cs_real, etpa_cs_imag, etpa_prob, tpa_cs)

    def calc_etpa_dw(self, nstate_in='all', wp=3.0996, ist=0, fst_in='auto', erange=1.5, kap=0.0, Te=4000, Ae=1e-06, tstep=1, line='l', pol=[-1, 4, -1], dw=0.1):
        if self.nst_calc:
            nstate = nstate_in
            fst = fst_in
        else:
            nstate, fst = self.calc_nst(self, nstate_in, wp, fst_in, erange)
        F = pol[0]
        G = pol[1]
        H = pol[2]
        w0 = wp / 2.0
        etpa_cs_tot = []
        etpa_cs_real = []
        etpa_cs_imag = []
        tpa_cs = 0.0
        etpa_prob = []

        dw_step = 0.01
        dw_list = [-dw + dw_step*i for i in range(0,int((2*dw+dw_step)/dw_step))]

        delta_tpa = 1 / (0.1 * math.pi)
        S_tpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                 [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                 [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
        tau1 = Te * 1e-15
        tau2 = Te * 1e-15
        denom = 1.0 / (30 * Ae * Te * 1e-15)

        for i, d in enumerate(dw_list):
            etpa_cs_tot.append(0.0)
            etpa_cs_real.append(0.0)
            etpa_cs_imag.append(0.0)
            etpa_prob.append(0.0)
            S_etpa = [[complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                      [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)],
                      [complex(0.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]]
            for f in fst:
                delta = cnst.ev2hz / self.calc_em_tau(0,f)
                for n in range(0, nstate):
                    delE1 = abs(self.states[n].energy - (w0 + d))
                    delE2 = abs(self.states[n].energy - (w0 - d))
                    if delE1 == 0:
                        continue
                    #    delE1 = float(1e-30)
                    #if delE2 == 0:
                    #    delE2 = float(1e-30)
                    if f == n:
                        tdip1 = [ self.interst_tdip[f][n][k] - self.interst_tdip[ist][ist][k] for k in range(0, 3) ]
                    else:
                        tdip1 = self.interst_tdip[f][n]
                    tdip2 = self.interst_tdip[n][ist]
                    for a in range(0, 3):
                        for b in range(0, 3):
                            S_etpa[a][b] += tdip1[a] * tdip2[b] * (1 - math.e ** ((complex(0.0, -1.0) * tau1 * delE1 - 0.5 * tau1 * kap) * cnst.ev2hz)) / (delE1 - complex(0.0, 0.5) * kap) + tdip1[b] * tdip2[a] * (1 - math.e ** ((complex(0.0, -1.0) * tau2 * delE2 - 0.5 * tau2 * kap) * cnst.ev2hz)) / (delE2 - complex(0.0, 0.5) * kap)
                            if i == 0:
                                S_tpa[a][b] += tdip1[a] * tdip2[b] / (delE1 - complex(0.0, 0.5) * kap) + tdip1[b] * tdip2[a] / (delE2 - complex(0.0, 0.5) * kap)

            for a in range(0, 3):
                for b in range(0, 3):
                    etpa_cs_real[i] += (F * S_etpa[a][a].real * S_etpa[b][b].real + G * S_etpa[a][b].real * S_etpa[a][b].real + H * S_etpa[a][b].real * S_etpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_cs_imag[i] += (F * S_etpa[a][a].imag * S_etpa[b][b].imag + G * S_etpa[a][b].imag * S_etpa[a][b].imag + H * S_etpa[a][b].imag * S_etpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[i] += (F * S_tpa[a][a].real * S_tpa[b][b].real + G * S_tpa[a][b].real * S_tpa[a][b].real + H * S_tpa[a][b].real * S_tpa[b][a].real) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    etpa_prob[i] += (F * S_tpa[a][a].imag * S_tpa[b][b].imag + G * S_tpa[a][b].imag * S_tpa[a][b].imag + H * S_tpa[a][b].imag * S_tpa[b][a].imag) * cnst.etpa_fact * w0 ** 2 * delta * denom
                    if i == 0:
                        tpa_cs += (2 * S_tpa[a][a].real * S_tpa[b][b].real + 2 * S_tpa[a][b].real * S_tpa[a][b].real + 2 * S_tpa[a][b].real * S_tpa[b][a].real) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30
                        tpa_cs += (2 * S_tpa[a][a].imag * S_tpa[b][b].imag + 2 * S_tpa[a][b].imag * S_tpa[a][b].imag + 2 * S_tpa[a][b].imag * S_tpa[b][a].imag) * cnst.tpa_fact * w0 ** 2 * delta_tpa / 30

            etpa_cs_tot[i] = etpa_cs_real[i]

        if len(etpa_cs_tot) >= 0:
            self.tpa_cs_calc = True
            self.etpa_cs_calc = True
        else:
            print 'Error: ETPA cross section not calculated'
        return (
         dw_list, etpa_cs_tot, etpa_cs_real, etpa_cs_imag, etpa_prob, tpa_cs)

    def write_etpa_cs(self, etpa_matrix, outfilename):
        Te = etpa_matrix[0]
        etpa_cs_tot = etpa_matrix[1]
        etpa_cs_real = etpa_matrix[2]
        etpa_cs_imag = etpa_matrix[3]
        etpa_prob = etpa_matrix[4]
        tpa_cs = etpa_matrix[5]
        etpaout = open(outfilename + '.etpa', 'w')
        etpaout.write('# Classical TPACS (cm^4s-1) = %s\n' % string.ljust(('{:.3e}').format(tpa_cs), 7))
        for t in range(0, len(Te)):
            etpaout.write('%s  %s  %s  %s  %s\n' % (str(Te[t]).rjust(5), string.rjust(('{:.3e}').format(etpa_cs_tot[t]), 7),
             string.rjust(('{:.3e}').format(etpa_cs_real[t]), 7), string.rjust(('{:.3e}').format(etpa_cs_imag[t]), 7),
             string.rjust(('{:.3e}').format(etpa_prob[t]), 7)))

        etpaout.close()

    def exc_lrtz(self, state_n='all', emax=10.0, emin=0.0, gamma=0.2, prog='mopac', infile=''):
        if not self.states_read:
            self.read_states()

        nstate = 0
        espacing = 0.005

        while self.states[(nstate + 1)].energy < emax:
            nstate += 1
            if nstate == len(self.states) - 1:
                break
	
        nstate = nstate if state_n == 'all' else min(state_n, nstate)
        self.calc_interst_tdip(nstate + 1, prog, infile, True)
        absint = [0.0] * int((emax - emin) / espacing + 1)
        stick_out = open('exc_stick.txt', 'w')
        lrtz_out = open('exc_lrtz.txt', 'w')
        stick_out.write('# E/eV  E/nm     f\n')
        lrtz_out.write('# E/eV     E/nm     abs/a.u.\n')
        for i in range(1, nstate + 1):
            stick_out.write('%s %s %s\n' % (string.rjust('%.2f' % self.states[i].energy, 6), string.rjust('%.2f' % float(1239.84187 / self.states[i].energy), 7),
             string.rjust('%.3f' % self.interst_tdip[0][i][3], 7)))
            for j in range(0, len(absint)):
                dE = j * espacing + emin - self.states[i].energy
                absint[j] += self.interst_tdip[0][i][3] * gamma / ((dE ** 2 + gamma ** 2) * math.pi)

        for i in range(0, len(absint)):
            E = (i + 0.5) * espacing + emin
            lrtz_out.write('%s %s %s\n' % (string.rjust('%.2f' % E, 6), string.rjust('%.2f' % float(1239.84187 / E), 7),
             string.rjust('%.3f' % absint[i], 7)))

        stick_out.close()
        lrtz_out.close()

    def singlet_lrtz(self, state_n='all', emax=5.0, emin=0.0, gamma=0.2, prog='mopac', infile=''):
        if not self.states_read:
            self.read_states()
        nstate = 0
        espacing = 0.01
        while self.states[(nstate + 1)].energy < emax:
            nstate += 1
            if nstate == len(self.states) - 1:
                break

        nstate = nstate if state_n == 'all' else min(state_n, nstate)
        self.calc_interst_tdip(nstate + 1, prog, infile)
        absint = [
         0.0] * int((emax - emin) / espacing + 1)
        stick_out = open('singlet_stick.txt', 'w')
        lrtz_out = open('singlet_lrtz.txt', 'w')
        stick_out.write('# E/eV  E/nm     f\n')
        lrtz_out.write('# E/eV     E/nm     abs/a.u.\n')
        for si in range(0, nstate + 1):
            for sj in range(si + 1, nstate + 1):
                Eij = self.states[sj].energy - self.states[si].energy
                if Eij != 0.0:
                    stick_out.write('%s %s %s\n' % (string.rjust('%.2f' % Eij, 6), string.rjust('%.2f' % float(1239.84187 / Eij), 7),
                     string.rjust('%.3f' % self.interst_tdip[si][sj][3], 7)))
                    for i in range(0, len(absint)):
                        dE = i * espacing + emin - Eij
                        absint[i] += self.interst_tdip[si][sj][3] * gamma / ((dE ** 2 + gamma ** 2) * math.pi)

        for i in range(0, len(absint)):
            E = (i + 0.5) * espacing + emin
            lrtz_out.write('%s %s %s\n' % (string.rjust('%.2f' % E, 6), string.rjust('%.2f' % float(1239.84187 / E), 7),
             string.rjust('%.3f' % absint[i], 7)))

        print 'Stick spectrum: singlet_stick.txt'
        print 'Lorentzian spectrum: singlet_lrtz.txt'
        stick_out.close()
        lrtz_out.close()

    def calc_em_tau(self, istate=0, fstate=1, prog='nwchem', infile=''):
        if not self.states_read:
            self.read_states()
        if not self.interst_tdip_calc:
            self.calc_interst_tdip(fstate + 1, prog, infile, onlyground=False)
        tau_inv = 0.0
        Ef = self.states[fstate].energy
        for i in range(istate, fstate):
            Ei = self.states[i].energy
            dE = abs(Ef - Ei) * cnst.ev2hz
            dip = 0.0
            for k in range(0, 3):
                dip += (self.interst_tdip[i][fstate][k] * cnst.au2ang * 1e-10) ** 2

            tau_inv += 4.0 * cnst.alpha * dE ** 3 * dip / (3.0 * cnst.speed ** 2)

        return tau_inv

    def get_mo_fock(self, init=0, final=1):
        if not self.orbs_read:
            self.read_orbs()
        init = self.homo - 2
        final = self.homo + 2
        self.mo_fock = []
        for i in range(init, final):
            self.mo_fock.append([])
            for j in range(init, final):
                elem = self.mos[i].energy if i == j else 0.0
                self.mo_fock[(i - init)].append(elem)
                print elem,

            print

    def rot_matrix(self, init=0, final=1, ind1=0, ind2=1, gam=0.0):
        U_dim = final - init
        U = []
        for i in range(U_dim):
            U.append([])
            for j in range(U_dim):
                elem = 0.0
                if i == ind1 and j == ind1:
                    elem = math.cos(gam)
                else:
                    if i == ind1 and j == ind2:
                        elem = math.sin(gam)
                    else:
                        if i == ind2 and j == ind1:
                            elem = -math.sin(gam)
                        else:
                            if i == ind2 and j == ind2:
                                elem = math.cos(gam)
                U[i].append(elem)

        return U

    def rotate_mo(self, mo1, mo2, gam):
        new_mo1 = mo1
        new_mo2 = mo2
        for i in range(len(mo1.coeff)):
            new_mo1.coeff[i] = math.cos(gam) * mo1.coeff[i] + math.sin(gam) * mol2.coeff[i]
            new_mo2.coeff[i] = -math.sin(gam) * mo1.coeff[i] + math.cos(gam) * mol2.coeff[i]

        return (
         new_mo1, new_mo2)
