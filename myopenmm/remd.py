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
inpf = '../inpdir/' + stage + '/nvt' + index + '.inp'
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_0 = md_nvt.setup(sysgro, systop)

# EM
simulation_em = md_nvt.minimize(simulation_0, 'em', index,  mddir='MD/')

# NVT
simulation_nvt = md_nvt.mdrun(simulation_em, 'nvt', index, mddir='MD/')

# REMD (NPT)
T_list = [453, 463, 473, 483]#, 493, 503, 513, 523 ]
n_replica = len(T_list)
remd = REMDConductor(T_list)
inpf = '../inpdir/' + stage + '/npt' + index + '.inp'
remd.loadFile(inpf)
nvtgro, systop = remd.preparation(sysdir='SYS/', mddir='MD/')
simulation_npt = remd.setup(nvtgro, systop)

# REMD
arglists = remd.setSims(simulation_npt, 'npt', index)
sims, enes = [], []

for i, args in enumerate(arglists):
    print('index:', i)
    j = '{0:02d}'.format(i)
    sim, energy = remd.remdrun(args[0], args[1], args[2], args[3], args[4])
    sims.append(sim)
    enes.append(energy)
print('enes:', enes)

