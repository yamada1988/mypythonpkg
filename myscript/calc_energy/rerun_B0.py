from incropenmm import *
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
gpuindex = args[3]

md_npt = MDConductor(index)
inpf = '../inpdir/' + stage + '/npt_B0.inp'

# Production Run
md_npt.loadFile(inpf)
sysgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD')
simulation_npt = md_npt.setup(sysgro, systop, deviceindex=gpuindex)

context = simulation_npt.context


ofname = 'DAT/npt' + index + '_B0.dat'
fname = 'MD/npt' + index + '.xtc'
#fname = 'MD/energy' + index + '.trr'
with open(ofname, 'wt') as f:
    f.write('# step\tEnergy(kJ/mol)\n')

# load trajcetory in chunk_traj
i = 1
for chunk_traj in md.iterload(fname,chunk=500,top=sysgro):
    # read box informations
    ucell_vs = chunk_traj.unitcell_vectors

    for j in range(len(chunk_traj)):
        # read and set positions
        context.setPositions(chunk_traj.openmm_positions(j))
        # set periodicbixvectors
        ucell_v = ucell_vs[j]

        context.setPeriodicBoxVectors(ucell_v[0],ucell_v[1],ucell_v[2])
        state = context.getState(getEnergy=True)
        energyval = state.getPotentialEnergy()
        energyval /= kilojoules/mole
        #print(i, energyval)
        with open(ofname, 'a+') as f:
            f.write('{0:04d}\t{1:10.3f}\n'.format(i, energyval))    
        i += 1
