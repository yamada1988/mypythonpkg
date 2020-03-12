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
numiter = args[4]

nframe =     10
niter  =    100
recstep =    10
totaltime = nframe * niter * recstep

md_nvt = MDConductor(index)
inpf = '../inpdir/' + stage + '/nvt_ghost.inp'
traj = 'MD/slt' + index + '.xtc'

# Original
md_nvt.loadFile(inpf)
sysgro, systop = md_nvt.preparation(sysdir='SYS/', mddir='MD')
simulation_nvt = md_nvt.setup(sysgro, systop, deviceindex=gpuindex)
context = simulation_nvt.context
positions = context.getState(getPositions=True).getPositions()

# Refs
simsys = context.getSystem()
numatm = simsys.getNumParticles()
masses = [simsys.getParticleMass(i) for i in range(numatm)]
totalMass = sum(masses)
for i in range(numatm):
    if i not in md_nvt.ghost_particles:
        simsys.setParticleMass(i, 0)


# Set velocities
def set_vel(context_):
    velocities = context_.getState(getVelocities=True).getVelocities()
    #print(velocities[0])
    for i in range(numatm):
        if i not in md_nvt.ghost_particles:
            velocities[i] = Vec3(0.0, 0.0, 0.0) * nanometer/picosecond
    context_.setVelocities(velocities)
    #print([simsys.getParticleMass(0)])

# Recreate simulation
topology = simulation_nvt.topology
integrator = VerletIntegrator(md_nvt.dt*picosecond)
simulation_refs = Simulation(topology, simsys, integrator, md_nvt.pltform, md_nvt.properties)
context = simulation_refs.context
context.setPositions(positions)
context.setVelocitiesToTemperature(md_nvt.temperature)

# XTC
mdxtc = 'MD/refs' + index + '_' + numiter + '.xtc'
xtc_reporter = XTCReporter(mdxtc, recstep)
simulation_refs.reporters.append(xtc_reporter)

# Log
mdlog = 'MD/refs' + index + '_' + numiter + '.log'
log_f = open(mdlog, mode='a+')
log_reporter = StateDataReporter(log_f, recstep, time=True,
                                                 totalEnergy=True, potentialEnergy=True, temperature=True, density=True,
                                                 progress=True, remainingTime=True, speed=True, totalSteps=totaltime, systemMass=totalMass, separator='\t')
simulation_refs.reporters.append(log_reporter)

# Run
it0 = (int(numiter) - 1) * nframe
itN = it0 + nframe
ts = md.load(traj, top=sysgro)
ucell_vs = ts.unitcell_vectors
for it,t in enumerate(ts):
    if it < it0:
        continue
    elif it == itN:
        break
    # set velocities
    set_vel(context)
    # set positions
    pos = t.openmm_positions(0)
    print('################')
    print('frame = {0:04d}'.format(it))
    print(pos[0])
    context.setPositions(pos)
    # set periodicbixvectors
    ucell_v = ucell_vs[it]
    context.setPeriodicBoxVectors(ucell_v[0],ucell_v[1],ucell_v[2])
    simulation_refs.step(recstep*niter) 
    #print([simsys.getParticleMass(0)])
    velocities = context.getState(getVelocities=True).getVelocities()
    print(velocities[0])
