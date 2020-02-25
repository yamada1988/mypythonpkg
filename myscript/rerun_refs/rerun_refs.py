from myopenmm import *
import mdtraj as md
import parmed as parm
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


nframe = 5000
niter  = 10
recstep = 20
totaltime = nframe * niter * recstep

md_nvt = MDConductor(index)
inpf = '../inpdir/' + stage + '/nvt_ghost.inp'
traj = 'MD/md' + index + '.xtc'

# Original
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_nvt = md_nvt.setup(sysgro, systop, deviceindex=gpuindex)
context = simulation_nvt.context

# Refs
simsys = simulation_nvt.system
numatm = simsys.getNumParticles()
masses = [simsys.getParticleMass(i) for i in range(numatm)]
totalMass = sum(masses)
for i in range(numatm):
    if i not in md_nvt.ghost_particles:
        simsys.setParticleMass(i, 0.0E0)

# XTC
mdxtc = 'MD/refs' + index + '.xtc'
xtc_reporter = XTCReporter(mdxtc, recstep)
simulation_nvt.reporters.append(xtc_reporter)


# Log
mdlog = 'MD/refs' + index + '.log'
log_f = open(mdlog, mode='a+')
log_reporter = StateDataReporter(log_f, recstep, time=True,
                                                 totalEnergy=True, potentialEnergy=True, temperature=True, density=True,
                                                 progress=True, remainingTime=True, speed=True, totalSteps=totaltime, systemMass=totalMass, separator='\t')
simulation_nvt.reporters.append(log_reporter)

# Run
ts = md.load(traj, top=sysgro)
ucell_vs = ts.unitcell_vectors
for it,t in enumerate(ts):
    # set positions
    pos = t.openmm_positions(0)
    context.setPositions(pos)
    # set periodicbixvectors
    ucell_v = ucell_vs[it]
    context.setPeriodicBoxVectors(ucell_v[0],ucell_v[1],ucell_v[2])
    simulation_nvt.step(recstep*niter) 
