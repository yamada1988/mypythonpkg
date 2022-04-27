from myopenmm import *
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
inpf = '../inpdir/' + stage + '/nvt_ghost' + index + '.inp'
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_0 = md_nvt.setup(sysgro, systop)

# EM
simulation_em = md_nvt.minimize(simulation_0, 'em', index,  sysdir='SYS/')

# NVT
nvtlog = 'MD/nvt' + index + '.log'
simulation_nvt, _ = md_nvt.mdrun(simulation_em, 'nvt', index, mddir='MD/')


inpf = '../inpdir/' + stage + '/npt_ghost' + index + '.inp'
niter = 1
# REST (NPT)

N=6
T_list = ['topol_{0:02d}.top'.format(i) for i in range(1,N+1)]
remd = REMDConductor(T_list, mode='REST')
# equilibriation
remd.conduct(inpf, index, equilibriation=True, parallel=False)

# sampling
inpf2 = '../inpdir/stage2' + '/npt_ghost' + index + '.inp'
niter = 1000
remd.conduct(inpf2, index, niter, equilibriation=False, parallel=True)
