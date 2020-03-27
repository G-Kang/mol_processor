#!/usr/bin/python

import sys,os,glob
from util import *
from numpy import *

freq = raw_input("Frequency as written in output name:")
mdir = "mode"+freq+"m"
pdir = "mode"+freq+"p"

os.system("tar -zxvf {0}.tar.gz".format(mdir))
os.system("tar -zxvf {0}.tar.gz".format(pdir))

data = open('gridinfo','r')

grid = []
for line in data.readlines():
    if '#' in line.split()[0]: continue
    grid.append([float(i) for i in line.split()])

gridn = []
for k in range(0,3):
    gridn.append((grid[k][1]-grid[k][0])/grid[k][2] + 1)

x = linspace(grid[0][0],grid[0][1],gridn[0])
y = linspace(grid[1][0],grid[1][1],gridn[1])
z = linspace(grid[2][0],grid[2][1],gridn[2])

currdir = os.getcwd()

for base in [mdir,pdir]:
    os.chdir(currdir+"/"+base)

    inname = base+'_grid.in'
    adfin = open(inname,'w')

    adfin.write("$ADFBIN/densf -n $NSLOTS << eor\n\nGrid Inline\n")

    for i in x:
	for j in y:
	    for k in z:
		adfin.write("%.2f\t%.2f\t%.2f\n" % (i,j,k))

    adfin.write("End\nDensity scf\naoresponse alpha lifetime\neor\n\n")
    adfin.close()
    os.system("~/bin/adf.sh {0}".format(inname))
