import math

test

#fosc_fact = 0.0038
#alph = 0.000029376024  # Used in polarizability calculation

# Energy units
ev2cm  = 8065.54
au2ev  = 27.2107
au2cm  = 219474.63
au2ang = 0.5291772
ev2j   = 1.60218E-19
fs2ev  = 1.519267515    # fs to eV-1
ev2hz  = 1.0/4.135667516*1.0E15
nm2ev  = 1.0/1240

ang2m = 1.0E-10
au2s  = 6.582119E-16

planck = 6.6260755E-34
hbar = planck/(2*math.pi)
boltz = 1.3806580E-23
speed = 2.99792458E+8
epsilon0 = 8.8541878E-12
avogadro = 6.02214199E+23
ea2debye = 4.80321
debye2cm = 1E-08/ea2debye
debye2au = 0.393430307
au2debye = 2.541746
ea2au = ea2debye*debye2au
me = 9.10938356E-31
amu = 1/avogadro*1E-03
temper = 300.0
exparg = planck*speed*100.0/boltz
conver = 2*planck*math.pi**2/speed*1.0E-40/amu
width = 20.0
scale = 1E+34
alpha = 1/137.035999679

afac = width/(2.0*math.pi)
bfac = width/2.0

## Natural Unit Conversion
nu2cm = hbar/ev2j*speed*100
au2nu = (au2ang*1E-8)/nu2cm
ev2nu = ev2j/hbar
fs2nu = ev2hz*1E-15
tpa_fact = (math.pi/2)*(nu2cm)**4
etpa_fact = (math.pi/4)*(nu2cm)**2

##
fosc_fact = 2.0/(3.0*au2ev)
#tpa_fact = 8*math.pi**3*(au2ang*1.0E-8)**4*alpha**2
#etpa_fact = 4*math.pi**3*(au2ang*1.0E-8)**2*alpha**2

# Polarizability unit conversion
alph = au2cm/au2ang**2/ea2debye**2

# Atomic masses
mass_dict = {
    'H'  :   1.008, 
    'He' :   4.0026,
    'Li' :   6.94,
    'Be' :   9.0122,
    'B'  :  10.81,
    'C'  :  12.011,
    'N'  :  14.007,
    'O'  :  15.999,
    'F'  :  18.998,
    'Na' :  22.989,
    'Mg' :  24.305,
    'Al' :  26.981,
    'Si' :  28.085,
    'P'  :  30.973,
    'S'  :  32.06,
    'Cl' :  35.45,
    'K'  :  39.0983,
    'Ca' :  40.078,
    'Fe' :  55.845,
    'Co' :  58.933,
    'Ni' :  58.6934,
    'Cu' :  63.546,
    'Zn' :  65.38,
    'As' :  74.9216,
    'Se' :  78.96,
    'Br' :  79.904,
    'Ru' : 101.07,
    'Rh' : 102.90,
    'Pd' : 106.42,
    'Ag' : 107.8682,
    'Cd' : 112.411,
    'In' : 114.818,
    'Sn' : 118.71,
    'I'  : 126.90,
    'Ir' : 192.217,
    'Pt' : 195.084,
    'Au' : 196.96,
    'Hg' : 200.59
}


valelec_dict = {
    'H'  :  1, 
    'He' :  0,
    'Li' :  1,
    'Be' :  2,
    'B'  :  3,
    'C'  :  4,
    'N'  :  5,
    'O'  :  6,
    'F'  :  7,
    'Na' :  1,
    'Mg' :  2,
    'Al' :  3,
    'Si' :  4,
    'P'  :  5,
    'S'  :  6,
    'Cl' :  7,
    'K'  :  1,
    'Ca' :  2,
    'Zn' :  2,
    'Ag' : 11,
    'Au' : 11
}
