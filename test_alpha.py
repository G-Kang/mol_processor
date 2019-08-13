#!/usr/bin/python

import string, sys
import constants as cnst

# Constants for later
gamma = 0.1088j
energy = 1.0
tdip = 5.0
step = 0.01
max = 0.0
nsteps = int(max/step) + 1

for i in range(0,nsteps):
    omega = i*step

    # Compute contribution to alpha
    alpha = complex(0.0,0.0)

    a_num = tdip**2 
    # Following Reimers code, convert to cm-1 here
    a_denom_1 = (energy - gamma - omega) * cnst.ev2cm
    a_denom_2 = (energy + gamma + omega) * cnst.ev2cm

    alpha += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph/3
    print string.rjust('%.3f'%omega,9), string.rjust('%.4f'%alpha.real,9),string.rjust('%.4f'%alpha.imag,9), string.rjust('%.4f'%abs(alpha),9)

