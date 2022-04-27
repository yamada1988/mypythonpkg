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

md_nvt = MDConductor()
inpf = '../inpdir/' + stage + '/nvt' + index + '.inp'
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_0 = md_nvt.setup(sysgro, systop)

# EM
simulation_em = md_nvt.minimize(simulation_0, 'em', index,  mddir='MD/')

# NVT
nvtlog = 'MD/nvt' + index + '.log'
simulation_nvt, _ = md_nvt.mdrun(simulation_em, 'nvt', index, mddir='MD/')

inpf = '../inpdir/' + stage + '/npt' + index + '.inp'
niter = 1
# REMD (NPT)
T_list = [453, 456, 459, 461, 463]
remd = REMDConductor(T_list, mode='REMD')
# equilibriation
remd.conduct(inpf, index, equilibriation=True)

# sampling
inpf2 = '../inpdir/stage2' + '/npt' + index + '.inp'
niter = 100
remd.conduct(inpf2, index, niter, equilibriation=False)
