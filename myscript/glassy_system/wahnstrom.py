# example illustrating use of Free Energy OpenMM plugin.
# Calculation of chemical potential of argon in a periodic box.
#
# AUTHOR
#
# John D. Chodera <jchodera@berkeley.edu>
#
# REQUIREMENTS
#
# numpy - Scientific computing package - http://numpy.scipy.org
#
# REFERENCES
#
# [1] Michael R. Shirts and John D. Chodera. Statistically optimal analysis of samples from multiple equilibrium states.
# J. Chem. Phys. 129:124105 (2008)  http://dx.doi.org/10.1063/1.2978177

from __future__ import print_function
from simtk.openmm import *
from simtk.unit import *
import numpy
import parmed as pmd
import sys


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
            serial = i+1
            name = resName = atoms[i]
            chainID = 'A'
            resSeq = i+1
            iCode = ''
            occupancy = 1.0
            tempFactor = 0.0
            x, y, z = coordinates[i]
            xyz.write("ATOM  {0:5d} {1:4s} {2:3s} {3:1s}{4:4d}{5:1s}   {6:8.3f}{7:8.3f}{8:8.3f}{9:6.2f}{10:6.2f}\n".format(serial, name, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor))
            
            
print(""" ###########################
   WAHNSTROM BINARY SYSTEM 
 ###########################""")

# =============================================================================
# Specify simulation parameters
# =============================================================================

nparticles = 216 # number of total particles
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
nequil_steps = 50000 # number of dynamics steps for equilibration
nsample_steps = 50000 # number of dynamics steps for sampling

# temperature
args = sys.argv
temp = float(args[1])

# Argon units
kB = BOLTZMANN_CONSTANT_kB
temp0 = epsilon1 / (kB*AVOGADRO_CONSTANT_NA)
NA0 = AVOGADRO_CONSTANT_NA / mole
m0 = mass1
tau0 = ((m0/amu*1.661*10**(-27))*(sigma1/angstrom*10**-10)**2/(epsilon1/kilojoule/mole/NA0)/1000.0e0)**0.50 * 10**12 * picosecond
print(temp0, epsilon1, tau0)
print("Target temperature:")
Tstar = float(temp)
print(Tstar)

reduced_density = 0.775 # reduced density rho*sigma^3
temperature = Tstar * kelvin # temperature
collision_rate = 2.0 / picosecond # collision rate for Langevin thermostat
timestep = 2.0 * femtosecond # integrator timestep

# =============================================================================
# Compute box size.
# =============================================================================
volume = nparticles*(sigma1**3)/reduced_density
box_edge = volume**(1.0/3.0)
cutoff = min(box_edge*0.49, 4.0*sigma1) # Compute cutoff
print("sigma1, sigma2, box_edge, cutoff")
print(sigma1, sigma2, box_edge, cutoff)

# =============================================================================
# Build systems
# =============================================================================

# Create random initial positions.
import numpy.random
positions = Quantity(numpy.random.uniform(high=box_edge/angstroms, size=[nparticles,3]), angstrom)
pbcbox = [box_edge, box_edge, box_edge]
pbcbox = list(map(lambda x: x/angstrom, pbcbox))
pos = list(map(lambda x: x/angstrom, positions))
pdbf = 'SYS/system0001.pdb'
write_pdb(pdbf, atoms, pos, pbcbox)
pdb = pmd.load_file(pdbf)

# Create argon system where first particle is alchemically modified by lambda_value.
system = System()
system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))

# Retrieve the NonbondedForce
nbforce = NonbondedForce()
for particle_index in range(nparticles):
    if particle_index < n1:
        atoms[particle_index] = 'Ar'
        system.addParticle(mass1)
        # Add normal particle.
        nbforce.addParticle(charge1, sigma1, epsilon1)
    elif n1 <= particle_index  < nparticles:
        atoms[particle_index] = 'Br'
        system.addParticle(mass2)
        # Add normal particle.
        nbforce.addParticle(charge2, sigma2, epsilon2)
system.addForce(nbforce)

#for particle_index in range(nparticles):
#    c_, s_, e_ = nbforce.getParticleParameters(particle_index)
#    print(particle_index, c_, s_, e_)

# =============================================================================
# Run simulations at each alchemical lambda value.
# Reduced potentials of sampled configurations are computed for all alchemical states for use with MBAR analysis.
# =============================================================================

# Create Integrator and Context.
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.reporters.append(mdtraj.reporters.HDF5Reporter(mdh5, 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, time=True, step=True,
    potentialEnergy=True, temperature=True, density=True,progress=True, remainingTime=True,
    speed=True, totalSteps=totaltime, separator='\t'))

# Initiate from last set of positions.
simulation.context.setPositions(positions)

# Minimize energy from coordinates.
print("minimizing..." )
simulation.LocalEnergyMinimizer.minimize()

# Equilibrate.
print("equilibrating...")
simulation.step(nequil_steps)


# append reporters
totaltime = 2500000
mdh5 = 'MD/md0001_{0:03d}.h5'.format(int(temp))
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
simulation = Simulation(top.topology, system, integrator0, platform)
simulation.reporters.append(mdtraj.reporters.HDF5Reporter(mdh5, 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, time=True, step=True,
    potentialEnergy=True, temperature=True, density=True,progress=True, remainingTime=True,
    speed=True, totalSteps=totaltime, separator='\t'))

# run 500 ps of production simulation
print('Running Production...')
simulation.step(totaltime)

# output grofile
mdgro = 'SYS/out0001_{0:03d}.gro'.format(int(temp))
with open(sysgro, 'rt') as f:
    lines = [line.strip('\n') for line in f]
    anum = lines[1]
    lines_ = lines[2:]

positions = simulation.context.getState(getPositions=True).getPositions()
box = simulation.context.getState().getPeriodicBoxVectors()
pbcbox = [box[0][0], box[1][1], box[2][2]]
pos = [pos for pos in positions]
pbcbox = list(map(lambda x: x/nanometer, pbcbox))
pos = list(map(lambda x: x/nanometer, pos))
dt_now = datetime.datetime.now()
date_ = "-{0:%Y-%m-%d-%H}".format(dt_now)
with open(mdgro, 'wt') as f:
    f.write('Generated by OpenMM: Date = {0:%Y-%m-%d %H:%M:%S}\n'.format(dt_now))
    f.write(' '+anum+'\n')
    for i,line in enumerate(lines_[:-1]):
        l = line[:20] + '{0:8.4f}{1:8.4f}{2:8.4f}\n'.format(pos[i][0], pos[i][1], pos[i][2])
        f.write(l)
    f.write('{0:7.4f}\t{0:7.4f}\t{0:7.4f}\n'.format(pbcbox[0], pbcbox[1], pbcbox[2]))

print('Done!')
