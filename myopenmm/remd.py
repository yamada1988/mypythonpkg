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

#md_nvt = MDConductor()
#inpf = '../inpdir/' + stage + '/nvt' + index + '.inp'
#md_nvt.loadFile(inpf)
#sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
#simulation_0 = md_nvt.setup(sysgro, systop)

# EM
#simulation_em = md_nvt.minimize(simulation_0, 'em', index,  mddir='MD/')

# NVT
#nvtlog = 'MD/nvt' + index + '.log'
#simulation_nvt, _ = md_nvt.mdrun(simulation_em, 'nvt', index, mddir='MD/')

inpf = '../inpdir/' + stage + '/npt' + index + '.inp'
# REMD (NPT)
T_list = [453, 455, 457]
remd = REMDConductor(T_list)
remd.loadFile(inpf)
nvtgro, systop = remd.preparation(sysdir='SYS/', mddir='MD/')
simulation_npt = remd.setup(nvtgro, systop)
print('pbcbox:', simulation_npt.system.getDefaultPeriodicBoxVectors())
remd.mdrun(simulation_npt, 'npt', index)
simulations = remd.spread_replicas(simulation_npt)
tstates = remd.initialize_replicas(simulations)
# REMD
# equilibriation and initialization
tstates, energys = remd.remdrun(tstates, 'npt', index, mddir='MD/', sysdir='SYS/')
print('enes:', energys)
b_list, p_list = remd.calc_prob(energys, 1)
print('blist:', b_list, 'plist:', p_list)
print(tstates)
tstates = remd.exchange(tstates, b_list, p_list, 1, index)
print(tstates)
sys.exit()

for niter in range(1, 2):    

    b_list, p_list = remd.calc_prob(enes, niter)
    print('blist:', b_list, 'plist:', p_list)
    remd.exchange(sims, b_list, p_list, niter, index)
