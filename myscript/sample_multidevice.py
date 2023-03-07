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
ens      = args[3]

md = MDConductor(index)
inpf = '../../inpdir_multidevice/' + stage + '/' + ens + '.inp'

# Production Run
md.loadFile(inpf)
sysgro, systop = md.preparation(sysdir='SYS/', mddir='MD')
simulation = md.setup(sysgro, systop, deviceindex="0,1,2,3")

simulation, _ = md.mdrun(simulation, ens, index, mddir='MD/')
