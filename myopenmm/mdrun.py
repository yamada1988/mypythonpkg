import myopenmm as mymm
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

inpf = 'openmm.inp'
md = mymm.MDConductor()
md.loadFile(inpf)
md.conduct(sysdir='SYS/', emname='em', nvtname='nvt', nptname='npt', mdname='md')
