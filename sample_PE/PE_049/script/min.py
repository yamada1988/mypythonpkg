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

md_min = MDConductor(index)
inpf = '../inpdir/' + stage + '/min.inp'

# Production Run
md_min.loadFile(inpf)
sysgro, systop = md_min.preparation(sysdir='SYS/', mddir='MD')
simulation = md_min.setup(sysgro, systop, deviceindex=gpuindex)
simulation = md_min.minimize(simulation, 'min', index, 'SYS/', 'MD/')

