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
import cupy as cp
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
nparticles =   50000 
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
nequil_steps  =  5000 # number of dynamics steps for equilibration
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
    if force_max <= 1000.0e0:
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
atmstep =   5
t_ratio = fhstep/atmstep
logfile = 'MD/log0001_{0:03d}.txt'.format(int(temp))
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
platform = Platform.getPlatformByName('OpenCL')
properties = {'Precision':'single'}
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
          'dt_fh/dt_atm': t_ratio}
dt_atom = atmstep*timestep * 10**-3 /picosecond

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
    lm = np.zeros((rN,rN,rN))
    lg = np.zeros((3,rN,rN,rN))
    ls0k = np.zeros((3,3,rN,rN,rN))
    lsv = np.zeros((3,3,rN,rN,rN))
    local_mass = 0.0e0
    local_g = 0.0e0
    local_v = 0.0e0
    local_sigma = 0.0e0
    local_sigma0k = 0.0e0
    local_sigmav = 0.0e0
    for t in range(fhstep/atmstep):
        simulation.step(atmstep)
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
        if t == 0:
            r_ = np.array([r_min+dr*ir for ir in range(rN)])
            rax = np.array([r_, r_, r_])
        # calculate rho, g, v, sigma_ab, tau_ab, S_ab, eta as local variables
        w = cp.asarray(masses, dtype=np.float32)
        lm = cp.zeros((rN,rN,rN), dtype=cp.float32)
        r0 = cp.asarray(pos[:,0]/dr, dtype=cp.int32)
        r1 = cp.asarray(pos[:,1]/dr, dtype=cp.int32)
        r2 = cp.asarray(pos[:,2]/dr, dtype=cp.int32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, float32 w',
            'raw float32 lm',
            '''
              int ind[] = {r0,r1,r2};
              atomicAdd(&lm[ind], w);
            '''            
            )(r0,r1,r2,w,lm)
        lm = cp.asnumpy(lm)
        lm /= dV
        local_mass += lm
        del lm

        w = cp.asarray(masses*vel.T, dtype=np.float32)
        lg = cp.zeros((3,rN,rN,rN), dtype=cp.float32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, raw float32 w',
            'raw float32 lg',
            '''
            for (int j = 1; j < 4; j++) {
            int ind[] = {j,r0,r1,r2};
            int ind2[] = {j,i};
            atomicAdd(&lg[ind], w[ind2]);
            }
            '''
            )(r0,r1,r2,w,lg)
        lg = cp.asnumpy(lg)

        lg /= dV
        local_g += lg
        del lg

        dum = masses*vel.T
        s = dum.T[:,np.newaxis,:]*vel[:,:,np.newaxis]
        w = cp.asarray(s, dtype=cp.float32)
        ls0k = cp.zeros((3,3,rN,rN,rN), dtype=cp.float32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, raw float32 w',
            'raw float32 ls0k',
            '''
            for (int j = 1; j < 4; j++) {
              for (int k = 1; k < 4; k++){
                int ind[] = {j,k,r0,r1,r2};
                int ind2[] = {i,j,k};
                atomicAdd(&ls0k[ind], w[ind2]);
              }
            }
            '''
            )(r0,r1,r2,w,ls0k)
        ls0k = cp.asnumpy(ls0k) /dV
        local_sigma0k += ls0k
        del ls0k   

        dum = pos[:,np.newaxis,:]*frc[:,:,np.newaxis]
        w = cp.asarray(dum, dtype=cp.float32)
        lsv = cp.zeros((3,3,rN,rN,rN), dtype=cp.float32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, raw float32 w',
            'raw float32 lsv',
            '''
            for (int j = 1; j < 4; j++) {
              for (int k = 1; k < 4; k++){
                int ind[] = {j,k,r0,r1,r2};
                int ind2[] = {i,j,k};
                atomicAdd(&lsv[ind], w[ind2]);
              }
            }
            '''
            )(r0,r1,r2,w,lsv)
        lsv = cp.asnumpy(lsv) /dV
        local_sigmav += lsv
        del lsv 

    # Time average
    local_mass /= t_ratio
    local_g /= t_ratio
    local_sigma /= t_ratio
    local_sigma0k /= t_ratio
    local_sigmav /= t_ratio
    local_sigma0 = local_sigma0k + local_sigmav

    local_mass = cp.asnumpy(local_mass)
    local_g = cp.asnumpy(local_g)
    local_sigma0 = cp.asnumpy(local_sigma0)
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
