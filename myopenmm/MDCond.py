from myopenmm import *
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os
from joblib import Parallel, delayed

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

B_factor = []
#for ene in enes:
#    B_ factor.append() 
nu-G02[15:19:19]:~/calspa/OpenMM/POLYMER_MIXTURE/PMMA/PMMA_050/MELT_050>cat ~/software/mypythonpkg/myopenmm/MDCond.py
import myopenmm as mymm
import mdtraj as md
from distutils.util import strtobool
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import parmed as pmd
from sys import stdout
import os
import os.path
import sys 
from mypool import MyPool

#  Requirement:
#  python 2.7, openmm, mdtraj, parmed
#

def unwrap_remdrun(arg, **kwarg):
    return REMDConductor.remdrun(*arg, **kwarg)

class MDConductor:
    def __init__(self):
        self.InpDict = {'integrator': 'Langevin', 
                        'temperature': 300.0e0, 
                        'pressure': 1.0e0, 
                        'sysfile' : 'system.gro',
                        'forcefield': 'topol.top', 
                        'fileformat': None, 
                        'friction_const': 1, 
                        'dt': 0.0020, 
                        'mdensemble' : 'npt',
                        'setPME' : False,
                        'pme-alpha' : 0.026, 
                        'pme-nx' : 64,
                        'pme-ny' : 64,
                        'pme-nz' : 64,
                        'pme-ftol' : 1.0e-5,     
                        'nvtstep' : 0, 
                        'nptstep' : 0, 
                        'mdstep' :  0, 
                        'nvteqrecstep' : 0, 
                        'npteqrecstep' : 0,
                        'mdrecstep' : 0, 
                        'emflag' : False, 
                        'nvtflag' : False, 
                        'nptflag' : False, 
                        'nvtrecflag' : False, 
                        'nptrecflag' : False, 
                        'mdrecflag' :  False, 
                        'simlogflag' : True, 
                        'simlogstep' : 500,
                        'verbose' : True, 
                        'forcerecflag' : False,
                        'remdflag' : False, 
                        'pbc' : True,
                        'annealingflag' : False,
                        'reporterformat' : None, 
                        'constraints' : None, 
                        'nonbonded_method' : 'PME', 
                        'nonbonded_switchflag' : True,
                        'nonbonded_switch' : 1.0e0, 
                        'nonbonded_cutoffflag' : True,
                        'nonbonded_cutoff' : 1.20e0, 
                        'heavyhydrogen' : False, 
                        'removeCMMotion' : False, 
                        'ghost_particle' : False,
                        'path_ndxfile' : None,
                        'platform' : 'CUDA', 
                        'precision' : 'double',
                        'recflag': False}
        
        print('Initialize...')               
    def loadFile(self, file):
        try:
            f = open(file, 'rt')
        except IOError:
            sys.exit('Input file {0} not founded.'.format(file))

        with open(file, 'rt') as f:
            InpDict = {line.split('=')[0].strip() : line.split('=')[1].strip() \
                           for line in f if "#" not in line and "=" in line}

        for k in InpDict:
            v = InpDict[k]
            if v == 'False':
                InpDict[k] = False
            elif v == 'True':
                InpDict[k] = True

        self.InpDict.update(InpDict)
        InpDict = self.InpDict 

        self.integrator = InpDict['integrator']
        self.temperature = float(InpDict['temperature'])
        self.pressure = float(InpDict['pressure'])
        self.sysfile = InpDict['sysfile']
        self.forcefield = InpDict['forcefield']
        self.fileformat = InpDict['fileformat']
        self.friction_const = float(InpDict['friction_const'])
        self.dt = float(InpDict['dt'])
        self.mdensemble = InpDict['mdensemble']
        self.nvtstep = int(InpDict['nvtstep'])
        self.nptstep = int(InpDict['nptstep'])
        self.mdstep = int(InpDict['mdstep'])
        self.setPME = InpDict['setPME']
        self.pme_alpha = float(InpDict['pme-alpha'])
        self.pme_nx = int(InpDict['pme-nx'])
        self.pme_ny = int(InpDict['pme-ny'])
        self.pme_nz = int(InpDict['pme-nz'])
        self.pme_ftol = float(InpDict['pme-ftol'])
        self.nvteqrecstep = int(InpDict['nvteqrecstep'])
        self.npteqrecstep = int(InpDict['npteqrecstep'])
        self.mdrecstep = int(InpDict['mdrecstep'])
        self.emflag = InpDict['emflag']
        self.nvtflag = InpDict['nvtflag']
        self.nptflag = InpDict['nptflag']
        self.nvtrecflag = InpDict['nvtrecflag']
        self.nptrecflag = InpDict['nptrecflag']
        self.mdrecflag = InpDict['mdrecflag']
        self.simlogflag = InpDict['simlogflag']
        self.simlogstep = int(InpDict['simlogstep'])
        self.verbos = InpDict['verbose']
        self.forcerecflag = InpDict['forcerecflag']
        self.pbc = InpDict['pbc']
        self.remdflag = InpDict['remdflag']
        self.annealingflag = InpDict['annealingflag']
        self.switchflag = InpDict['nonbonded_switchflag']
        self.nonbonded_switch = float(InpDict['nonbonded_switch'])
        self.reporterformat = InpDict['reporterformat']
        self.constraints = InpDict['constraints']
        self.nonbonded_method = InpDict['nonbonded_method']
        self.nonbonded_cutoffflag = InpDict['nonbonded_cutoffflag']
        self.nonbonded_cutoff = float(InpDict['nonbonded_cutoff'])
        self.heavyhydrogen = InpDict['heavyhydrogen']
        self.removeCMMotion = InpDict['removeCMMotion'] 
        self.ghost_particle = InpDict['ghost_particle']
        self.path_ndxfile = InpDict['path_ndxfile']
        self.platform = InpDict['platform']
        self.precision = InpDict['precision']
        self.recflag = InpDict['recflag']

        with open('openmm.out', 'wt') as f:
            print('writing openmm.out...')
            for k, v in InpDict.items():
                k = k.ljust(20)
                f.write('{0}:{1}\n'.format(k, v))


    def preparation(self, sysdir='SYS/', mddir='MD/'):
        
        print('preparation start...')
        # Platform input
        pltfmname = self.platform
        precision = self.precision
        self.pltform = Platform.getPlatformByName(pltfmname)
        if self.platform == 'CUDA':
            self.properties = {'CudaPrecision': precision}#, 'CudaDeviceIndex': '0,1'}#, 'CudaUseBlockingSync':False }
        elif self.platform == 'CPU':
            self.properties = {}
        elif self.platform == 'OpenCL':
            self.properties = {'OpenCLPrecision': precision, 'OpenCLDeviceIndex': '0,1,2,3'}


        # System input
        sysgro = sysdir + self.sysfile
        systop = sysdir + self.forcefield

        # System output
        mddir = 'MD/'

        return sysgro, systop

    def setup(self, sysgro, systop):
        # Simulation setting
        if self.temperature:
            temperature = self.temperature * kelvin

        if self.friction_const:
            fric_const = self.friction_const / picosecond

        dt = self.dt * picosecond

        if self.pressure:
            pressure = self.pressure * bar

        # Create gro
        root, ext = os.path.splitext(sysgro)
        if ext == '.gro':
            gro = GromacsGroFile(sysgro)
        elif ext == '.pdb':
            gro = PDBFile(sysgro) 

        # Create top
        if self.pbc:
            print('set periodic boundary condition...')
            if ext == '.gro':
                top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            if ext == 'gro':
                top = GromacsTopFile(systop)
       
        if self.heavyhydrogen:
            print('Hydrogen Mass Repartitioning...')
            hmass = 4*amu
        else:
            hmass = 1*amu


        # Create system
        if self.nonbonded_method == 'PME':
            print('set PME...')
            if self.nonbonded_cutoffflag:
                nonbonded_cutoff = self.nonbonded_cutoff   
                system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=PME, 
                                          nonbondedCutoff=nonbonded_cutoff,
                                          constraints=self.constraints)
            else:
                system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=PME,
                                          constraints=self.constraints)        
        elif self.nonbonded_method == 'Cutoff':
            print('set cutoff...')
            nonbonded_cutoff = self.nonbonded_cutoff
            system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=Cutoff, 
                                      nonbondedCutoff=nonbonded_cutoff,
                                      constraints=self.constraints)

        # Check ensembleflag
        if self.nvtflag:
            nvtstep = self.nvtstep  
        else:
            nvtstep = 0

        if self.nptflag:
            nptstep = self.nptstep
            print('add MonteCarlo Barostat...')
            system.addForce(MonteCarloBarostat(pressure, temperature))
        else:
            nptstep = 0       

        self.steps = nvtstep + nptstep 

        if self.nvtrecflag:
            nvteqrecstep = self.nvteqrecstep
            self.recflag = True
        else:
            nvteqrecstep = 0
        if self.nptrecflag:
            npteqrecstep = self.npteqrecstep
            self.recflag = True
        else:
            npteqrecstep = 0

        if self.mdrecflag:
            mdrecstep =  self.mdrecstep
            self.recflag = True
        else:
            mdrecstep = 0
        self.recstep = nvteqrecstep + npteqrecstep + mdrecstep

        if self.simlogflag:
            simlogstep = self.simlogstep


        if self.setPME:
            print('set PME parameters...')
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            forces['NonbondedForce'].setPMEParameters(self.pme_alpha, self.pme_nx, self.pme_ny, self.pme_nz)
            print('after:')
            alpha, nx, ny, nz = forces['NonbondedForce'].getPMEParameters()
            print(alpha)
            print('nx={0:d}, ny={1:d}, nz={2:d}'.format(nx, ny, nz))

        if self.switchflag:
            print('set switching function...')
            nonbonded_switch = self.nonbonded_switch
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            forces['NonbondedForce'].setUseSwitchingFunction(True)
            forces['NonbondedForce'].setSwitchingDistance(nonbonded_switch)

        # Create integrator
        if self.integrator == 'Langevin':
            print('set Langevin Integrator...')
            integrator = LangevinIntegrator(temperature, fric_const, dt) 
        elif self.integrator == 'Brownian':
            print('set Brownian Integrator...')
            integrator = BrownianIntegrator(temperature, fric_const, dt)
        elif self.integrator == 'Verlet':
            print('set Verlet Integrator...')
            integrator = VerletIntegrator(dt)
        else:
            sys.exit('Invalid Integrator type. Check your input file.')
        print('dt:', dt, 'fric_const:', fric_const, 'Total step:', self.steps)
        
        # Ghost-particle        
        if self.ghost_particle:
            print('set ghost particle calculation...')
            inpf = self.path_ndxfile        
            with open(inpf, 'rt') as gf:
                total_lines =[line.strip() for line in gf]
                try:
                    core_index = total_lines.index('[ Core ]') + 1
                    ghost_index = total_lines.index('[ Ghost ]') + 1
                    solvent_index = total_lines.index('[ Solvent ]') + 1
                except:
                    sys.exit('If ghost_particle flag is True, ghost_index file must be specified.\n')
            core_start = total_lines[core_index].split('-')[0]
            core_end = total_lines[core_index].split('-')[1]
            ghost_start = total_lines[ghost_index].split('-')[0]
            ghost_end = total_lines[ghost_index].split('-')[1]
            solvent_start = total_lines[solvent_index].split('-')[0]
            solvent_end = total_lines[solvent_index].split('-')[1]
 
            if int(core_start) == 0 and int(core_start) == int(core_end):
                core_particles = range(0,0)
            else:
                core_particles = range(int(core_start)-1, int(core_end))
            ghost_particles = range(int(ghost_start)-1, int(ghost_end))
            solvent_particles = range(int(solvent_start)-1, int(solvent_end))

            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            nonbonded_force = forces['NonbondedForce']

            solute_particles = core_particles + ghost_particles
            slt_param = [[]]*len(solute_particles) 
            for index in solute_particles:
                slt_param[index] = nonbonded_force.getParticleParameters(index)
                if index in ghost_particles:
                    nonbonded_force.setParticleParameters(index, charge=0.0e0, sigma=slt_param[index][1], epsilon=0.0e0)

            for i in ghost_particles:
                for j in solute_particles:
                    if i == j:
                        continue
                    try:
                        nonbonded_force.addException(i, j, chargeProd=slt_param[i][0]*slt_param[j][0], 
                                                           sigma=0.50*(slt_param[i][1]+slt_param[j][1]),
                                                           epsilon=sqrt(slt_param[i][2]*slt_param[j][2]))
                    except:
                        pass
#                        print('{0:d}-{1:d} pair already prepared.'.format(i, j))
#            for index in range(nonbonded_force.getNumExceptions()):
#                print(nonbonded_force.getExceptionParameters(index))
            
        # Create simulation
        simulation = Simulation(top.topology, system, integrator, self.pltform, self.properties)
        simulation.context.setPositions(gro.positions)
        return simulation

    # EM simulation
    def minimize(self, simulation, emname, index, mddir='MD'):
        emname = emname + index

        print('Minimizing...')
        empdb = mddir + emname + '.pdb'
        simulation.minimizeEnergy(maxIterations=2000)

        print('Check Energy...')
        state = simulation.context.getState(getEnergy=True)
        energyval = state.getPotentialEnergy()
        print(energyval)
        eneval = energyval / (kilojoule/mole)
        if eneval > 1.0e+14:
            sys.exit('Energy diverges at infinity.')

        print('Saving...')
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open(empdb, 'w'))

        return simulation

    def mdrun(self, simulation, ensname, index, mddir='MD/', sysdir='SYS/'):
        # ensname simulation
        print('self.steps:', self.steps, 'mdresctep:', self.recstep)
        print(ensname+' Simulation...')
        mdname = ensname + index
        mdlog = mddir + mdname + '.log'
        simulation.reporters.append(StateDataReporter(mdlog, self.recstep, time=True,
                                                      totalEnergy=True, temperature=True, density=True, 
                                                      progress=True, remainingTime=True, speed=True, totalSteps=self.steps, separator='\t'))
        if self.recflag:
            print('\nSaving...')
            mdxtc = mddir + mdname + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, self.recstep)
            simulation.reporters.append(xtc_reporter)

        print('Check Energy...')
        state = simulation.context.getState(getEnergy=True)
        energyval = state.getPotentialEnergy()
        print(energyval)
        if energyval / (kilojoule/mole) > 1.0e+14:
            print('Energy diverged at infinity.')
            sys.exit()

        simulation.step(self.steps)
        print('Done!\n')

        mdpdb = sysdir + mdname + '.pdb'
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open(mdpdb, 'w'))

        root, ext = os.path.splitext(mdpdb)
        outf = root + '.gro'
        pdbstructure = pmd.load_file(mdpdb)
        pdbstructure.save(outf, format='GRO', overwrite=True)
        return simulation


    def convert_pdb2gro(self, pdb):
        fname = pdb
        root, ext = os.path.splitext(pdb)
        outf = root + '.gro'
        pdbstructure = pmd.load_file(pdb)
        pdbstructure.save(outf, format='GRO', overwrite=True)


    def getIntegratorInfo(self, integrator):
            dt = self.dt
            fric_const = integrator.getFriction()
            temperature = integrator.getTemperature()

            return dt, fric_const, temperature


# Below description mey only be valid for Python 2.X. 
# Modify ...(class, object) to ...(class) and super(class ,self) to super() if error occures.
class REMDConductor(MDConductor, object):
    def __init__(self, T_list):
        super(REMDConductor, self).__init__()
        self.T_list = list(map(float, T_list))

    def setSims(self, simulation, ensname, index, T_list, mddir='MD/', sysdir='SYS/'):
        T_list = self.T_list
        n_replica = len(T_list)
        print(n_replica, T_list)
        simulations = [ '' for i in range(n_replica)]

        dt, fric_const, temperature = self.getIntegratorInfo(simulation.integrator)
        print(dt, fric_const, temperature)
        system = simulation.system
        top = simulation.topology
        positions = simulation.context.getState(getPositions=True).getPositions()

        arglists = [ []*5 for line in range(n_replica)]
        for i, T_ in enumerate(T_list):
            integrator_ = LangevinIntegrator(T_, fric_const, dt)
            properties = {'Precision': self.precision}#, 'DeviceIndex': '{0:d}'.format(i)}
            simulations[i] = Simulation(top, system, integrator_, self.pltform,  properties)   
            simulations[i].context.setPositions(positions)
            j = '{0:02d}'.format(i+1)
            arglists[i] = [simulations[i], ensname, index+'_'+ j, mddir, sysdir]
        print(arglists)

        return arglists


    def wrapper(self, args):
        return self.mdrun(*args)


    def remdrun(self, simulation, ensname, index, mddir='MD/', sysdir='SYS/'):
        super(REMDConductor, self).mdrun(simulation, ensname, index, mddir='MD/', sysdir='SYS/')
        print('Check Energy for REMD...')
        state = simulation.context.getState(getEnergy=True)
        energyval = state.getPotentialEnergy()
        return simulation, energyval 
