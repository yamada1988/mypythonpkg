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

parm = pmd.load_file('SYS/topol_B0.top')

md_npt = MDConductor(index)
inpf = '../inpdir/' + stage + '/B0.inp'

# Production Run
md_npt.loadFile(inpf)
sysgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD')
simulation_npt = md_npt.setup(sysgro, systop)
#context = simulation_npt.context
#print(pmd.openmm.energy_decomposition(parm, context))



print("total:")
state = simulation_npt.context.getState(getEnergy=True)
energyval = state.getPotentialEnergy()
print(energyval)
