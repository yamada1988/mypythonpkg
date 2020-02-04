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
inpf = '../inpdir/' + stage + '/npt_ghost.inp'

# Production Run
md_npt.loadFile(inpf)
sysgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD')
simulation_npt = md_npt.setup(sysgro, systop)
context = simulation_npt.context


ofname = 'energy' + index + '_ghost.dat'
fname = 'MD/md' + index + '_01.xtc' 
with open(ofname, 'wt') as f:
    f.write('# step\tEnergy(kJ/mol)\n')

traj = md.load(fname,top=sysgro)
for i in range(len(traj)):
    context.setPositions(traj.openmm_positions(i))
    state = context.getState(getEnergy=True)
    energyval = state.getPotentialEnergy()
    energyval /= kilojoules/mole
    print(i, energyval)
    with open(ofname, 'a+') as f:
        f.write('{0:04d}\t{1:10.3f}\n'.format(i, energyval)) 
