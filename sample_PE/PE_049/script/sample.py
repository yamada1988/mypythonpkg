from incropenmm import *
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os

args = sys.argv
stage = args[1]
index = args[2]
gpuindex = args[3]


md_npt = MDConductor(index)
inpf = '../inpdir/' + stage + '/npt_ghost.inp'

# Production Run
md_npt.loadFile(inpf)
sysgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD')
simulation_npt = md_npt.setup(sysgro, systop, deviceindex=gpuindex)

simulation_npt, _ = md_npt.mdrun(simulation_npt, 'npt', index, mddir='MD/')
