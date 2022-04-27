# example illustrating use of on-the-fly calculation for local variables in OpenMM plugin.
# Calculation of Wahnstrom binary solution systems
#
# AUTHOR
#
# Kazuo Yamada
#
# REQUIREMENTS
#
# np - Scientific computing package - http://np.scipy.org
# h5py - Pythonic interface to the HDF5 binary data format - https://www.h5py.org
#
# REFERENCES
#
# [1] Michael R. Shirts and John D. Chodera. Statistically optimal analysis of samples from multiple equilibrium states.
# J. Chem. Phys. 129:124105 (2008)  http://dx.doi.org/10.1063/1.2978177
# [2] Goran Wahnstrom. Molecular-dynamics study of a supercooled two-component Lennard- Jones system.
# Phys. Rev. A. 44:6 (1991)

from __future__ import print_function
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import mdtraj
import numpy as np
import parmed as pmd
import h5py
import sys
from sys import stdout

def write_pdb(filename, atoms, coordinates, box, mode="w"):
    """Very primitive PDB writer.
    :Arguments:
       *filename*
           name of the output file
       *atoms*
           list of the N atom names
       *coordinates*
           coordinates as Nx3 array (must be in Angstroem)
       *box*
           box lengths (Lx Ly Lz) (must be in Angstroem)
    See http://www.wwpdb.org/documentation/format32/sect9.html
    """
    with open(filename, mode) as xyz:
        xyz.write("HEADER    simple PDB file with {0:d} atoms\n".format(len(atoms)))
        xyz.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}  90.00  90.00  90.00 P 1        1\n".format(box[0], box[1], box[2]))
        for i in xrange(len(atoms)):
            serial = (i+1) % 10000
            name = resName = atoms[i]
            chainID = 'A'
            resSeq = (i+1) % 10000
            iCode = ''
            occupancy = 1.0
            tempFactor = 0.0
            x, y, z = coordinates[i]
            xyz.write("ATOM  {0:5d} {1:4s} {2:3s} {3:1s}{4:4d}{5:1s}   {6:8.3f}{7:8.3f}{8:8.3f}{9:6.2f}{10:6.2f}\n".format(serial, name, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor))
            

# number of total particles            
nparticles = 500000 
print(""" ###########################
   WAHNSTROM BINARY SYSTEM 
   Num: {0:9d}
 ###########################""".format(nparticles))

# =============================================================================
# Specify simulation parameters
# =============================================================================

n1 = int(nparticles/2)
n2 = n1
atoms = ['' for i in range(nparticles)]
for particle_index in range(nparticles):
    if particle_index < n1:
        atoms[particle_index] = 'Ar'
    elif n1 <= particle_index  < nparticles:
        atoms[particle_index] = 'Br'

# species1
mass1 = 39.9 * amu # mass
sigma1 = 3.4 * angstrom # Lennard-Jones sigma
epsilon1 = 0.994 * kilojoule/mole # Lennard-Jones well-depth
charge1 = 0.0e0 * elementary_charge # argon model has no charge

# species2
mass2 = 79.896 * amu # mass
sigma2 = 4.08 * angstrom # Lennard-Jones sigma
epsilon2 = 0.994 * kilojoule/mole# Lennard-Jones well-depth
charge2 = 0.0e0 * elementary_charge # argon model has no charge

# Steps
nequil_steps  = 50000 # number of dynamics steps for equilibration
nsample_steps = 50000 # number of dynamics steps for sampling
rN = 32

# temperature
args = sys.argv
temp = float(args[1])

# Argon units
kB = BOLTZMANN_CONSTANT_kB
temp0 = epsilon1 / (kB*AVOGADRO_CONSTANT_NA)
NA0 = AVOGADRO_CONSTANT_NA / mole
m0 = mass1
l0 = sigma1
tau0 = ((m0/amu*1.661*10**(-27))*(sigma1/angstrom*10**-10)**2/(epsilon1/kilojoule/mole/NA0)/1000.0e0)**0.50 * 10**12 * picosecond
f0 = m0*l0/(tau0**2)
e0 = epsilon1
print(temp0, epsilon1, tau0)
Tstar = float(temp)
print("Target temperature:", Tstar)
kBT = kB * Tstar * kelvin * AVOGADRO_CONSTANT_NA / (kilojoule/mole) 

reduced_density = 0.750 # reduced density rho*sigma^3
temperature = Tstar * kelvin # temperature
collision_rate = 2.0 / picosecond # collision rate for Langevin thermostat
timestep = 2.0 * femtosecond # integrator timestep

# =============================================================================
# Compute box size.
# =============================================================================
volume = (n1*sigma1**3+n2*sigma2**3)/reduced_density
box_edge = volume**(1.0/3.0)
cutoff = min(box_edge*0.49, 2.50*sigma1) # Compute cutoff
print("sigma1, sigma2, box_edge, cutoff")
print(sigma1, sigma2, box_edge, cutoff)

# =============================================================================
# Build systems
# =============================================================================

# Create random initial positions
gridnum = int(box_edge/(2.0*angstrom))
print('total_divnum:',gridnum**3)
x_ = np.linspace(0, box_edge/angstrom, num=gridnum)
x_ = Quantity(x_ , angstrom)
x_, y_, z_ = np.meshgrid(x_, x_, x_)
xyzmesh = np.vstack((x_.flatten(), y_.flatten(), z_.flatten())).T

id_ = np.random.choice(gridnum**3,nparticles,replace=False)
xyz_ = np.round(np.array([xyzmesh[i]+1.0*np.random.rand(3) for i in id_]), 6)
positions = Quantity(xyz_ , angstrom)

pbcbox = [box_edge, box_edge, box_edge]
pbcbox = list(map(lambda x: x/angstrom, pbcbox))
pos = np.array(list(map(lambda x: x/angstrom, positions)))
pdbf = 'SYS/system0001.pdb'
write_pdb(pdbf, atoms, pos, pbcbox)
pdb_ = PDBFile(pdbf)

# Create argon system where first particle is alchemically modified by lambda_value.
system = System()
system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))

# Retrieve the NonbondedForce
print('Set nonbonded parameters...')
nbforce = NonbondedForce()
for particle_index in range(n1):
    atoms[particle_index] = 'Ar'
    system.addParticle(mass1)
    # Add normal particle.
    nbforce.addParticle(charge1, sigma1, epsilon1)
for particle_index in range(n1,n1+n2):
    atoms[particle_index] = 'Br'
    system.addParticle(mass2)
    # Add normal particle.
    nbforce.addParticle(charge2, sigma2, epsilon2)
nbforce.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
nbforce.setCutoffDistance(cutoff) 
system.addForce(nbforce)
print('Finish to set nonbonded parameters')

# Create Integrator and Context.
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
platform = Platform.getPlatformByName('CPU')
properties = {'Precision':'mixed'}
mdh5 = 'MD/mim0001_{0:03d}.h5'.format(int(temp))
mdh5_ = mdtraj.reporters.HDF5Reporter(mdh5, 1000)
simulation = Simulation(pdb_.topology, system, integrator, platform)
simulation.reporters.append(mdh5_)
simulation.reporters.append(StateDataReporter(stdout, 1000, time=True, step=True,
      potentialEnergy=True, temperature=True, density=True,progress=True, remainingTime=True,
      speed=True, totalSteps=nequil_steps, separator='\t'))

# Initiate from last set of positions.
simulation.context.setPositions(positions)

print('minimizing...')
for attmpt in range(1000):
    state = simulation.context.getState(getEnergy=True)
    energyval = state.getPotentialEnergy()
    print('attempt {0:6d}:'.format(attmpt))
    print('Enrg: {0:8.5e} (kJ/mol)'.format(energyval/(kilojoule/mole)))

    # Minimize energy from coordinates.
    simulation.minimizeEnergy(maxIterations=2)
    forces = simulation.context.getState(getForces=True).getForces(asNumpy=True) / (kilojoule/mole/nanometers) 
    force_max = np.amax(forces)
    print('Fmax: {0:8.5e} (kJ/mol/nm)'.format(force_max))
    if force_max <= 500.0e0:
        break
    else:
        mdh5_.close()

positions = simulation.context.getState(getPositions=True).getPositions()
box = simulation.context.getState().getPeriodicBoxVectors()
pbcbox = [box[0][0], box[1][1], box[2][2]]
pos = [pos for pos in positions]
pbcbox = list(map(lambda x: x/angstrom, pbcbox))
pos = list(map(lambda x: x/angstrom, pos))
write_pdb(pdbf, atoms, pos, pbcbox)
pdb_ = PDBFile(pdbf)
system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))
# Create Integrator and Context.
recstep =  50
fhstep  = 500
logfile = 'MD/log0001_{0:03d}.txt'.format(int(temp))
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
platform = Platform.getPlatformByName('OpenCL')
properties = {'Precision':'mixed'}
mdh5 = 'MD/md0001_{0:03d}.h5'.format(int(temp))
mdh5_ = mdtraj.reporters.HDF5Reporter(mdh5, recstep, potentialEnergy=False, kineticEnergy=False, temperature=False, velocities=True)
simulation = Simulation(pdb_.topology, system, integrator, platform, properties)
simulation.reporters.append(mdh5_)
simulation.reporters.append(StateDataReporter(logfile, recstep, time=True, step=True,
      potentialEnergy=True, temperature=True, density=True,progress=True, remainingTime=True,
      speed=True, totalSteps=nsample_steps+nequil_steps, separator='\t'))

# Initiate from last set of positions.
simulation.context.setPositions(positions)
# Equilibrate.
print("equilibrating...")
simulation.step(nequil_steps)

# run 500 ps of production simulation
print('On-the-fly calculation...')
units_ = {'length': 'nanometers',
          'mass': 'amu',
          'time': 'picosecond',
          'force': 'kilojoule/mole/nanometer',
          'energy': 'kilojoule/mole',
          'dt_fh/dt_atm': fhstep}
dt_atom = timestep * 10**-3 /picosecond

with h5py.File('MD/local0001_{0:03d}.h5'.format(int(temp)), 'w') as f:
    f.create_group('Units')
    for k,v in units_.items():
        f['Units'].create_dataset(k, data=v)
    f.create_dataset('dt_atom', data=dt_atom)
    f.create_dataset('localMass', shape=(nsample_steps/recstep, rN, rN, rN))
    f.create_dataset('localStressTensor', shape=(nsample_steps/recstep, 3, 3, rN, rN, rN))
    f.create_dataset('localMomentum', shape=(nsample_steps/recstep, 3, rN, rN, rN))
    f.create_dataset('MeshGrid', shape=(3, rN+1))
    f.create_dataset('time', shape=(nsample_steps/recstep,))

masses = np.array([system.getParticleMass(indx)/amu  for indx in range(nparticles)]) 
for it in range(nsample_steps/fhstep):
    print('it:',it)
    time_ = (it+1)*fhstep 
    lg = np.zeros((3,rN,rN,rN))
    lv = np.zeros((3,rN,rN,rN))
    ls0k = np.zeros((3,3,rN,rN,rN))
    lsv = np.zeros((3,3,rN,rN,rN))
    lS = np.zeros((3,3,rN,rN,rN))
    lvo = np.zeros((3,rN,rN,rN))
    local_mass = 0.0e0
    local_g = 0.0e0
    local_v = 0.0e0
    local_sigma = 0.0e0
    local_sigma0k = 0.0e0
    local_sigmav = 0.0e0
    for t in range(fhstep):
        print('t_atom:', t)
        simulation.step(1)
        pos = simulation.context.getState(getPositions=True).getPositions(asNumpy =True) / nanometers
        vel = simulation.context.getState(getVelocities=True).getVelocities(asNumpy =True) / (nanometers/picoseconds)
        frc = simulation.context.getState(getForces=True).getForces(asNumpy =True)   / (kilojoule/mole/nanometers)
        box = simulation.context.getState().getPeriodicBoxVectors(asNumpy =True) / nanometers
        L = box[0,0]
        # unset PBC 
        pos -= np.trunc(pos/L)*L
        pos += np.round((0.50e0*L-pos)/L)*L

        # grid parameters
        r_min = 0.0e0 
        r_max = L
        dr = (r_max-r_min)/ float(rN)
        dV = dr**3
        # print('dr:',dr, 'L:', L, 'dV:', dV)
        # calculate rho, g, v, sigma_ab, tau_ab, S_ab, eta as local variables
        lm, rax = np.histogramdd(pos, bins=(rN, rN, rN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max)), weights=masses)
        for xa in range(3):
            lg[xa], rax = np.histogramdd(pos, bins=(rN, rN, rN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max)), weights=masses*vel.T[xa])
        lm /= dV
        lg /= dV
        local_mass += lm
        local_g += lg
        for xa in range(3):
            for xb in range(3):
                sk = masses*vel.T[xa]*vel.T[xb]
                #print('mvv, kT:')
                #print(masses*vel.T[xa]*vel.T[xb], kBT)
                ls0k[xa][xb], rax = np.histogramdd(pos, bins=(rN, rN, rN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max)), weights=sk)
                #print('sigma_k:')
                #print(local_sigma_k[xa][xb])
                sv = pos.T[xa]*frc.T[xb]
                lsv[xa][xb], rax = np.histogramdd(pos, bins=(rN, rN, rN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max)), weights=sv)
                #print('sigma_v:')
                #print(local_sigma_v[xa][xb]/dV)
        lsv /= dV
        local_sigma0k += ls0k
        local_sigmav += lsv
   
    # Time average
    local_mass /= recstep
    local_g /= recstep
    local_sigma /= recstep
    local_sigma0k /= recstep
    local_sigmav /= recstep
    local_sigma0 = local_sigma0k + local_sigmav
    #print(local_sigmak, local_sigmav)

    meshgrid = np.array(rax)
    with h5py.File('MD/local0001_{0:03d}.h5'.format(int(temp)), 'a') as f:
        f['localMass'][it] = local_mass 
        f['localStressTensor'][it] = local_sigma0
        f['localMomentum'][it] = local_g
        if t == 0:
            f['MeshGrid'] = meshgrid 
        f['time'][it] = time_*timestep *10**-3 / femtosecond

# output pdbfile
mdpdb = 'SYS/out0001_{0:03d}.pdb'.format(int(temp))
positions = simulation.context.getState(getPositions=True).getPositions()
box = simulation.context.getState().getPeriodicBoxVectors()
pbcbox = [box[0][0], box[1][1], box[2][2]]
pos = [pos for pos in positions]
pbcbox = list(map(lambda x: x/angstrom, pbcbox))
pos = list(map(lambda x: x/angstrom, pos))
write_pdb(mdpdb, atoms, pos, pbcbox)

print('Done!')
