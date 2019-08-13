#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute Raman cross-sections from manually entered excited state info
#
# Rebecca Gieseking, 2/3/2017
# 
#-------------------------------------------------------------------------

import string,sys,math
import constants as cnst
from molecule import *
from vibrations import *
from properties import *
import getopt

###############
# Manual versions of scripts
###############
# compute_alpha
def comp_alpha(states, omega, gamma):
    alpha = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]

    for k in range(0,len(states)):
        s = states[k]
        # Compute contribution to alpha
        #a_state = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]
        for i in range(0,3):
            for j in range(0,3):
                a_num = s.tdip[i]*s.tdip[j]
                # Following Reimers code, convert to cm-1 here
                a_denom_1 = (s.energy - gamma - omega) * cnst.ev2cm
                a_denom_2 = (s.energy - gamma + omega) * cnst.ev2cm

                alpha[i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
                #a_state[i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph
        #a_or = (a_state[0][0] + a_state[1][1] + a_state[2][2]) / 3.0
    return alpha


# alpha_slope
def a_slope(mst, pst, disp, omega, gamma):
    ds = disp/cnst.bohr2ang

    alpha_minus = comp_alpha(mst, omega, gamma)
    alpha_plus  = comp_alpha(pst, omega, gamma)

    alpha_diff = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]
    for i in range(0,3):
        for j in range(0,3):
            alpha_diff[i][j] = (alpha_plus[i][j] - alpha_minus[i][j])*cnst.bohr2ang**2/(2.0*ds)
    return alpha_diff

# Compute the Raman scattering factor S
def ram_scat(coord, alpha_diff):
    if 'x' in coord or 'y' in coord or 'z' in coord:
        if 'x' in coord:
            j = 0
        elif 'y' in coord:
            j = 1
        elif 'z' in coord:
            j = 2
        asq_r = alpha_diff[j][j].real**2
        asq_i = alpha_diff[j][j].imag**2

        gsq_r = 0.0
        gsq_i = 0.0

        # Compare with complex number math
        asq_r = alpha_diff[j][j]**2
        gsq_r = complex(0.0,0.0)


    else:
        asq_r = ((alpha_diff[0][0].real + alpha_diff[1][1].real + alpha_diff[2][2].real)/3.0)**2
        asq_i = ((alpha_diff[0][0].imag + alpha_diff[1][1].imag + alpha_diff[2][2].imag)/3.0)**2

        gsq_r = (6.0 * (alpha_diff[0][1].real**2 + alpha_diff[1][2].real**2 + alpha_diff[2][0].real**2) \
              + (alpha_diff[0][0].real - alpha_diff[1][1].real)**2 \
              + (alpha_diff[1][1].real - alpha_diff[2][2].real)**2 \
              + (alpha_diff[2][2].real - alpha_diff[0][0].real)**2) / 2.0
        gsq_i = (6.0*(alpha_diff[0][1].imag**2 + alpha_diff[1][2].imag**2 + alpha_diff[2][0].imag**2) \
              + (alpha_diff[0][0].imag - alpha_diff[1][1].imag)**2 \
              + (alpha_diff[1][1].imag - alpha_diff[2][2].imag)**2 \
              + (alpha_diff[2][2].imag - alpha_diff[0][0].imag)**2) / 2.0

    s_fact_r = 45.0*asq_r + 7.0*gsq_r
    s_fact_i = 45.0*asq_i + 7.0*gsq_i

    #print string.rjust('%.5f'%(math.sqrt(s_fact_r**2 + s_fact_i**2)),14)

    return s_fact_r, s_fact_i 

# Compute the Raman differential cross-section
def ram_cross(omega, s_fact_r, s_fact_i, freq=1000):
    tempfc = 1.0e6/(1.0 - math.exp(-cnst.exparg*freq/cnst.temper))
    lambda0 = omega * cnst.ev2cm
    frq4th = (lambda0 - freq)**4

    crs_real = tempfc * frq4th * cnst.conver * s_fact_r / (45.0 * freq)
    crs_imag = tempfc * frq4th * cnst.conver * s_fact_i / (45.0 * freq)
    crs_tot = crs_real + crs_imag

    return crs_real, crs_imag, crs_tot

# Constants for later
types = []
omega = 0.0
gamma = 0.1088j
disp_str = '0.01'
disp = float(disp_str)
coord = ''
nstates = 1
width = 20.0
scale = 1E+33
outfilename = 'man_raman.out'

helpfile = """
manual_raman.py -s <# states> -o <output>

    -s    Number of states to include in SOS expression             Default = All states
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"ho:s:",['--help','--output=','--states='])
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
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-s','--states'):
        nstates = int(arg)

out = open(outfilename+'.out','w')

# Manually enter state info
out.write ("Raman calculations for manual excited states\n")
print "Enter information for %i states" % nstates

states = []
mstates = []
pstates = []
for i in range(0,nstates):
    print "State %i" % (i+1)
    e = float(raw_input("   Enter the excited state energy (eV):       "))
    o = float(raw_input("   Enter the oscillator strength:             "))
    p =       raw_input("   Enter the polarization (3 floats):         ")
    pol = p.split()
    for j in range(0,3):
        pol[j] = float(pol[j])
    s = float(raw_input("   Enter the energy change with displacement: "))

    states.append(State(e,o))
    states[i].calc_tdip(pol)

    mstates.append(State(e-s,o))
    mstates[i].calc_tdip(pol)

    pstates.append(State(e+s,o))
    pstates[i].calc_tdip(pol)

    out.write("State %i   Energy = %.2f eV    Osc = %.2f    Pol = %.2f %.2f %.2f    Slope = %.4f \n" % (i, e, o, pol[0], pol[1], pol[2], s))

# Manual version of Raman cross-sections
out.write("\n Energy (eV)               crs_r               crs_i               crs_t \n")
omega = 0.0
o_max = 5.0
o_step = 0.05
for i in range(0,int(o_max/o_step)+1):
    diff          = a_slope(mstates, pstates, disp, omega, gamma)
    s_r, s_i      = ram_scat(' ',diff)
    c_r, c_i, c_t = ram_cross(omega, s_r, s_i)

    out.write(string.rjust('%.4f'%omega,12) + string.rjust(str(c_r),20) + string.rjust(str(c_i),20) + string.rjust(str(c_t),20) + '\n')

    omega += o_step

out.close()


