from myopenmm import *
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os
import parmed as pmd

args = sys.argv
stage = args[1]
index = args[2]


md_npt = MDConductor(index)
inpf = '../inpdir/' + stage + '/S0.inp'

# Production Run
md_npt.loadFile(inpf)
sysgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD')
simulation_npt = md_npt.setup(sysgro, systop)
context = simulation_npt.context


ofname = 'energy' + index + '_openmm_S0.dat'
fname = 'MD/energy' + index + '.trr'

with open(ofname, 'wt') as f:
    f.write('# step\tEnergy(kJ/mol)\n')

# load trajcetory
traj = md.load(fname,top=sysgro)
# read box informations
ucell_vs = traj.unitcell_vectors

for i in range(len(traj)):
    # read and set positions
    context.setPositions(traj.openmm_positions(i))
    # set periodicbixvectors
    ucell_v = ucell_vs[i]

    context.setPeriodicBoxVectors(ucell_v[0],ucell_v[1],ucell_v[2])
    state = context.getState(getEnergy=True)
    energyval = state.getPotentialEnergy()
    energyval /= kilojoules/mole
    print(i, energyval)
    with open(ofname, 'a+') as f:
        f.write('{0:04d}\t{1:10.3f}\n'.format(i, energyval))    

