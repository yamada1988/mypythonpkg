import incropenmm as incrmm
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

md_nvt = incrmm.MDConductor()
inpf = '../inpdir/' + stage + '/nvt' + index + '.inp'
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_0 = md_nvt.setup(sysgro, systop)

# EM
simulation_em = md_nvt.minimize(simulation_0, 'em', index,  mddir='MD/')

# NVT
simulation_nvt, nvtpdb = md_nvt.mdrun(simulation_em, 'nvt', index, mddir='MD/')
md_nvt.convert_pdb2gro('SYS/nvt'+index+'.pdb')

# NPT
md_npt = incrmm.MDConductor()
inpf = '../inpdir/' + stage + '/npt' + index + '.inp'
md_npt.loadFile(inpf)
nvtgro, systop = md_npt.preparation(sysdir='SYS/', mddir='MD/')
simulation_npt = md_npt.setup(nvtgro, systop)
simulation_npt, nptpdb = md_npt.mdrun(simulation_npt, 'npt', index, mddir='MD/')
