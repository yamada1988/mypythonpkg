from __future__ import print_function
import mdtraj.reporters
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

print("Argon unit")
sigma0 = 0.34 * nanometer
ep0 = 0.99774 * kilojoule / mole
kB = BOLTZMANN_CONSTANT_kB
temp0 = ep0 / (kB*AVOGADRO_CONSTANT_NA)
NA0 = AVOGADRO_CONSTANT_NA / mole
m0 = 39.948 * 1.661e-27 * kilogram
tau0 = ((m0/kilogram)*(sigma0/nanometer*10**-9)**2/(ep0/kilojoule/mole/NA0)/1000.0e0)**0.50 * 10**12 * picosecond
print(temp0, ep0, tau0)
print("High temperature")
Tstar = temp0 *2.0e0
print(Tstar)

gro = app.GromacsGroFile('SYS/system0001.gro')
top = app.GromacsTopFile('SYS/topol.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=CutoffPeriodic, nonbondedCutoff=1*nanometer)
integrator = LangevinIntegrator(Tstar*kelvin, 1/picosecond, 0.002*picoseconds)

# prepare simulation
platform = Platform.getPlatformByName('OpenCL')
properties = {'Precision': 'mixed'}

#simulation = Simulation(top.topology, system, integrator, platform, properties)
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)

# minimize
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=500000)

# equilibrate for 100 steps
simulation.context.setVelocitiesToTemperature(Tstar*kelvin)
print('Equilibrating...')
simulation.step(100)

# append reporters
totaltime = 10000
simulation.reporters.append(mdtraj.reporters.HDF5Reporter('MD/test.h5', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, time=True, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=totaltime, separator='\t'))

# run 500 ps of production simulation
print('Running Production...')
simulation.step(totaltime)

print('add Barostat...')
system.addForce(MonteCarloBarostat(1.0*atmospheres, Tstar*kelvin, 25))
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)
simulation.context.setVelocitiesToTemperature(Tstar*kelvin)
totaltime = 250000
simulation.reporters.append(mdtraj.reporters.HDF5Reporter('MD/test.h5', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, time=True, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
     speed=True, totalSteps=totaltime, separator='\t'))

print('Running Production...')
simulation.step(totaltime)

print('Done!')
