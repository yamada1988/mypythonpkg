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
import numpy as np
from numpy.random import *
import math
import threading
import subprocess
import shutil
import datetime

#  Requirement:
#  python 2.7, openmm, mdtraj, parmed
#
def get_gpu_info():
   cmd = 'nvidia-smi --query-gpu=index --format=csv'
   output = subprocess.check_output(cmd, shell=True)
   lines = output.split('\n')
   lines = [ line.strip() for line in lines if line.strip() != '' ]
   print('GPU IDs:', lines)
   return lines

class MDRUN:
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
                        'nvestep' : 0,
                        'mdstep' :  0, 
                        'nverecstep': 0,
                        'nvteqrecstep' : 0, 
                        'npteqrecstep' : 0,
                        'mdrecstep' : 0, 
                        'emflag' : False, 
                        'nveflag': False,
                        'nvtflag' : False, 
                        'nptflag' : False, 
                        'nverecflag': False,
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
                        'recflag': False,
                        'virtualflag': False,
                        'multidevice': None,
                        'flatbottom': False,
                        'restraint_file': None}
        
        print('Initialize...')               
    def loadFile(self, file):
        try:
            f = open(file, 'rt')
        except IOError:
            sys.exit('Input file {0} not found.'.format(file))

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
        self.nvestep = int(InpDict['nvestep'])
        self.mdstep = int(InpDict['mdstep'])
        self.setPME = InpDict['setPME']
        self.pme_alpha = float(InpDict['pme-alpha'])
        self.pme_nx = int(InpDict['pme-nx'])
        self.pme_ny = int(InpDict['pme-ny'])
        self.pme_nz = int(InpDict['pme-nz'])
        self.pme_ftol = float(InpDict['pme-ftol'])
        self.nverecstep = int(InpDict['nverecstep'])
        self.nvteqrecstep = int(InpDict['nvteqrecstep'])
        self.npteqrecstep = int(InpDict['npteqrecstep'])
        self.mdrecstep = int(InpDict['mdrecstep'])
        self.emflag = InpDict['emflag']
        self.nveflag = InpDict['nveflag']
        self.nvtflag = InpDict['nvtflag']
        self.nptflag = InpDict['nptflag']
        self.nverecflag = InpDict['nverecflag']
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
        self.reporterformat = InpDict['reporterformat']
        self.virtualflag = InpDict['virtualflag']
        self.multidevice = InpDict['multidevice']
        self.flatbottom = InpDict['flatbottom']
        self.restraint_file = InpDict['restraint_file']

        with open('openmm.out', 'wt') as f:
            print('writing openmm.out...')
            for k, v in InpDict.items():
                k = k.ljust(20)
                f.write('{0}:{1}\n'.format(k, v))

    def create_platform(self, deviceindex=None):
        pltfmname = self.platform
        precision = self.precision
        self.pltform = Platform.getPlatformByName(pltfmname)
        if self.platform == 'CUDA':
            if deviceindex is None:
                if self.multidevice is None:
                    self.properties = {'Precision': precision}
                else:
                    self.properties = {'Precision': precision, 'DeviceIndex': self.multidevice}
            else:
                self.properties = {'Precision': precision, 'DeviceIndex': deviceindex}
        elif self.platform == 'CPU':
            self.properties = {}
        elif self.platform == 'OpenCL':
            if deviceindex is None:
                if self.multidevice is None:
                    self.properties = {'Precision': precision}
                else:
                    self.properties = {'Precision': precision, 'DeviceIndex': self.multidevice}
            else:
                self.properties = {'OpenCLPrecision': precision, 'DeviceIndex': deviceindex}
        print(self.properties)

    def preparation(self, sysdir='SYS/', mddir='MD/'):
        
        print('preparation start...')
        # System input
        sysgro = sysdir + self.sysfile
        systop = sysdir + self.forcefield

        with open(sysgro, 'rt') as f:
            lines = [line.strip('\n') for line in f]
            if self.fileformat == 'GRO':
                self.anum = lines[1]
                self.lines = lines[2:]
            elif self.fileformat == 'PDB':
                self.anum = [lines[1],lines[3]]
                self.lines = lines[4:]
        # System output
        mddir = 'MD/'

        print('sysgro:', sysgro, 'systop:', systop)
        return sysgro, systop


    def create_top(self, systop, gro):
        if self.pbc:
            print('set periodic boundary condition...')
            top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            top = GromacsTopFile(systop)

        self.top = top
        return top

    def create_system(self, top):
        # Create system
        if self.heavyhydrogen:
            print('Hydrogen Mass Repartitioning...')
            self.hmass = 4*amu
        else:
            self.hmass = None

        if self.nonbonded_method == 'PME':
            print('set PME...')
            if self.nonbonded_cutoffflag and self.fileformat == 'GRO':
                nonbonded_cutoff = self.nonbonded_cutoff   
                system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=PME, 
                                          nonbondedCutoff=self.nonbonded_cutoff,
                                          constraints=self.constraints,rigidWater=True)

            elif self.nonbonded_cutoffflag and self.fileformat == 'PDB':
                nonbonded_cutoff = self.nonbonded_cutoff   
                system = top.createSystem(self.modeller.topology, hydrogenMass=self.hmass,nonbondedMethod=PME, 
                                          nonbondedCutoff=self.nonbonded_cutoff,
                                          constraints=self.constraints,rigidWater=True)

            elif not self.nonbonded_cutoff_flag:
                system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=PME,
                                          constraints=self.constraints)        


        elif self.nonbonded_method == 'Cutoff':
            print('set cutoff...')
            system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=Cutoff, 
                                      nonbondedCutoff=self.nonbonded_cutoff,
                                      constraints=self.constraints)
        

        # Check ensembleflag
        if self.nveflag:
            nvestep = self.nvestep
        else:
            nvestep = 0

        if self.nvtflag:
            nvtstep = self.nvtstep  
        else:
            nvtstep = 0

        if self.pressure:
            pressure = self.pressure * bar

        if self.nptflag:
            nptstep = self.nptstep
            print('add MonteCarlo Barostat...')
            system.addForce(MonteCarloBarostat(pressure, self.temperature))
        else:
            nptstep = 0       

        self.steps = nvestep + nvtstep + nptstep 
        if self.nverecflag:
            nverecstep = self.nverecstep
            self.recflag = True
        else:
            nverecstep = 0
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
        self.recstep = nverecstep + nvteqrecstep + npteqrecstep + mdrecstep

        if self.simlogflag:
            simlogstep = self.simlogstep

        if self.setPME:
            print('set PME parameters...')
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            forces['NonbondedForce'].setPMEParameters(self.pme_alpha, self.pme_nx, self.pme_ny, self.pme_nz)
            #print('after:')
            alpha, nx, ny, nz = forces['NonbondedForce'].getPMEParameters()
            print(alpha)
            print('nx={0:d}, ny={1:d}, nz={2:d}'.format(nx, ny, nz))

        if self.switchflag:
            print('set switching function...')
            nonbonded_switch = self.nonbonded_switch
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            forces['NonbondedForce'].setUseSwitchingFunction(True)
            forces['NonbondedForce'].setSwitchingDistance(nonbonded_switch)

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

            solute_particles = list(core_particles) + list(ghost_particles)
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            for k, v in forces.items():
                print(k,v)
            nonbonded_force = forces['NonbondedForce']


            total_particles = system.getNumParticles()
            slt_param = [[]]*total_particles
            for index in range(total_particles):
                slt_param[index] = nonbonded_force.getParticleParameters(index)
                if index in ghost_particles:
                    nonbonded_force.setParticleParameters(index, charge=0.0e0, sigma=slt_param[index][1], epsilon=0.0e0)
            
            for i in core_particles:
                for j in core_particles:
                    if i == j:
                        continue
                    try:
                        nonbonded_force.addException(i, j,
                                                chargeProd = slt_param[i][0]*slt_param[j][0],
                                                sigma = (slt_param[i][1]+slt_param[j][1])/2.0E0,
                                                epsilon = np.sqrt(slt_param[i][2]*slt_param[j][2]))
                        #print(i, j)
                    except:
                        #print(i,j)
                        pass
            
            # NonbondedForce 2: Core-Ghost interactionm
            nonbonded_force_02 = NonbondedForce()
            
            cf = nonbonded_force.getCutoffDistance()
            me = nonbonded_force.getNonbondedMethod()
            tl = nonbonded_force.getEwaldErrorTolerance()
            sw = nonbonded_force.getSwitchingDistance()
            nonbonded_force_02.setCutoffDistance(cf)
            nonbonded_force_02.setNonbondedMethod(me)
            nonbonded_force_02.setPMEParameters(alpha, nx, ny, nz)
            nonbonded_force_02.setSwitchingDistance(sw)
            
            for index in range(total_particles):
                nonbonded_force_02.addParticle(1.0, 1.0, 1.0)
                if index in solute_particles:
                    nonbonded_force_02.setParticleParameters(index, charge=slt_param[index][0], sigma=slt_param[index][1], epsilon=slt_param[index][2])
                else:
                    nonbonded_force_02.setParticleParameters(index, charge=0.0e0, sigma=slt_param[index][1], epsilon=0.0e0)
            
            
            # Set Exculusion Info
            for index in range(nonbonded_force.getNumExceptions()):
                exception_info = nonbonded_force.getExceptionParameters(index)
                nonbonded_force_02.addException(exception_info[0], exception_info[1], 0.0e0, exception_info[3], 0.0e0)
                nonbonded_force_02.setExceptionParameters(index, exception_info[0], exception_info[1], 0.0e0, exception_info[3], 0.0e0)
            
            
           #system.addForce(nonbonded_force_02)



            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            #for k, v in forces.items():
            #     print(k, v)


        # Set Restraint
        # Refs: http://www.maccallumlab.org/news/2015/1/23/testing
        if self.flatbottom == True:   
            flat_bottom_force = CustomBondForce(
	    'step(r-r0) * (k/2) * (r-r0)^2')
            flat_bottom_force.addPerBondParameter('r0')
            flat_bottom_force.addPerBondParameter('k')

            with open(self.restraint_file) as input_file:
    	        for line in input_file:
                    columns = line.split()
                    atom_index_i = int(columns[0]) - 1
                    atom_index_j = int(columns[1]) - 1
                    k = float(columns[2])
                    r0 = float(columns[3])
            print('flat-bottom restraint added particle-{0:d} and particle-{1:d}.'.format(atom_index_i, atom_index_j))
            print('k:{0:8.5f}\tr0:{1:5.3f}'.format(k, r0))
            flat_bottom_force.addBond(
            atom_index_i, atom_index_j, [r0, k])

            system.addForce(flat_bottom_force)


        return system

   
    def create_integrator(self):
        if self.temperature:
            temperature = self.temperature * kelvin

        if self.friction_const:
            fric_const = self.friction_const / picosecond
        dt = self.dt * picosecond

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
        print('dt:', dt, 'fric_const:', fric_const, 'Temperature:', temperature, 'Total step:', self.steps)
        return integrator


    def setup(self, sysgro, systop, deviceindex=None,sysdir='SYS/'):
        # Create gro
        if self.fileformat == 'GRO':
            gro = GromacsGroFile(sysgro)
            # Create top
            top = self.create_top(systop, gro)

        elif self.fileformat == 'PDB':
            pdb = PDBFile(sysgro)
            # Create top
            top = ForceField(sysdir+self.forcefield)
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.addExtraParticles(top)
            self.modeller = modeller

        # Create system
        system = self.create_system(top)

        # Create integrator
        integrator = self.create_integrator()
       
        # Create platform
        self.create_platform(deviceindex) 
        print ('Create simulation...')
        if self.fileformat == 'GRO':
            simulation = Simulation(top.topology, system, integrator, self.pltform, self.properties)
            simulation.context.setPositions(gro.positions)
        elif self.fileformat == 'PDB':
            simulation = Simulation(self.modeller.topology, system, integrator, self.pltform, self.properties)
            print(simulation.topology)
            simulation.context.setPositions(self.modeller.positions)
        simulation.context.setVelocitiesToTemperature(self.temperature)

        return simulation

    # EM simulation
    def minimize(self, simulation, emname, index, sysdir='SYS/',max_iter=3500):
        emname = emname + index

        print('Minimizing...')
        empdb = sysdir + emname + '.pdb'
        simtop = simulation.topology
        simsys = simulation.system
        simint = LangevinIntegrator(300, 1.0, 2.0) 
        simpos = simulation.context.getState(getPositions=True).getPositions()
        simulation_ = Simulation(simtop, simsys, simint, Platform.getPlatformByName('CPU'))
        simulation_.context.setPositions(simpos)
        
        simulation_.minimizeEnergy(maxIterations=max_iter,tolerance = 4.148)

        print('Check Energy...')
        state = simulation_.context.getState(getEnergy=True)
        energyval = state.getPotentialEnergy()
        print(energyval)
        eneval = energyval / (kilojoule/mole)
        if eneval > 1.0e+14:
            sys.exit('Energy diverges at infinity.')

        print('Saving...')
        positions = simulation_.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation_.topology, positions, open(empdb, 'w'))

        simulation.context.setPositions(positions)
        return simulation

    def mdrun(self, simulation, ensname, index, mddir='MD/', sysdir='SYS/', 
              assert_system=True, check_eneflag=False, nrep='01', remdflag=False, niter=1, niter_tot=1):
        print('')
        print(ensname+' Simulation:', 'self.steps:', self.steps, 'mdresctep:', self.recstep)
        mdname = ensname + index
        if remdflag == False:
            mdlog = mddir + mdname + '.log'
        elif remdflag == True:
            index_new = index + '_' + nrep
            mdlog = mddir + ensname + index_new + '.log'
    
        log_f = open(mdlog, mode='a+')
        log_reporter = StateDataReporter(log_f, self.recstep, time=True,
                                                              totalEnergy=True, potentialEnergy=True, temperature=True, density=True,
                                                              progress=True, remainingTime=True, speed=True, totalSteps=niter_tot*self.steps, separator='\t')
        simulation.reporters.append(log_reporter)

        xtc_flag = False
        hdf5_flag = False
        if remdflag == True:
            print('nrep:', nrep, 'ID:', simulation.context.sysid)

        if self.recflag and self.reporterformat == 'XTC' and remdflag == False:
            xtc_flag = True
            #print('Save trajectory as xtcfile...')
            mdxtc = mddir + mdname + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, self.recstep)
            simulation.reporters.append(xtc_reporter)
        elif self.recflag and self.reporterformat == 'HDF5' and remdflag == False:
            hdf5_flag = True
            #print('Save trajectory as xtcfile...')
            mdh = mddir + mdname + '.h5'
            hdf5_reporter = md.reporters.HDF5Reporter(mdh, self.recstep, kineticEnergy=False,velocities=True)
            simulation.reporters.append(hdf5_reporter)
        elif self.recflag and simulation.context.sysid == '01':
            xtc_flag = True
            #print('Save only tagged trajectory in REMD simulation as xtcflie...')
            indexj = index + '_{0:04d}'.format(niter)
            mdname_new = ensname + indexj
            mdxtc = mddir + mdname_new + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, self.recstep)
            simulation.reporters.append(xtc_reporter)

        if assert_system == True:
            #print('Assert system...')
            state = simulation.context.getState(getEnergy=True)
            energyval = state.getPotentialEnergy()
            print(energyval)
            if energyval / (kilojoule/mole) > 1.0e+14:
                print('Energy diverged at infinity.')
                sys.exit()

        #print('Conducting...')
        simulation.step(self.steps)
        #print('Done!\n')


        # output grofile
        if niter == niter_tot:
            #print('Saving...')
            if remdflag == False:
                if self.fileformat == 'GRO':
                    mdgro = sysdir + mdname + '.gro'
                elif self.fileformat == 'PDB':
                    mdpdb = sysdir + mdname + '.pdb'
            elif remdflag == True:
                index_new = index + '_' + nrep
                mdgro = sysdir + ensname + index_new + '.gro'
            positions = simulation.context.getState(getPositions=True).getPositions()
            box = simulation.context.getState().getPeriodicBoxVectors()
            pbcbox = [box[0][0], box[1][1], box[2][2]]
            pos = [pos for pos in positions]
            pbcbox = list(map(lambda x: x/nanometer, pbcbox))
            pos = list(map(lambda x: x/nanometer, pos))
            dt_now = datetime.datetime.now()
            date_ = "-{0:%Y-%m-%d-%H}".format(dt_now)
            if self.fileformat == 'GRO':
                if os.path.exists(mdgro):
                    shutil.copyfile(mdgro, mdgro.split('.')[0] + date_ + '.gro')
                with open(mdgro, 'wt') as f:
                    f.write('OpenMM: Date = {0:%Y-%m-%d %H:%M:%S}\n'.format(dt_now))
                    f.write(' '+self.anum+'\n')
                    for i,line in enumerate(self.lines[:-1]):
                        l = line[:20] + '{0:8.4f}{1:8.4f}{2:8.4f}\n'.format(pos[i][0], pos[i][1], pos[i][2])
                        f.write(l)
                    f.write('{0:7.4f}\t{0:7.4f}\t{0:7.4f}\n'.format(pbcbox[0], pbcbox[1], pbcbox[2]))
            if self.fileformat == 'PDB':
                if os.path.exists(mdpdb):
                    shutil.copyfile(mdpdb, mdpdb.split('.')[0] + date_ + '.pdb')
                with open(mdpdb, 'wt') as f:
                    f.write('OpenMM: Date = {0:%Y-%m-%d %H:%M:%S}\n'.format(dt_now))
                    f.write(' '+self.anum[0]+'\n')
                    f.write('CRYST1  {0:7.3f}  {0:7.3f}  {0:7.3f}  90.00  90.00  90.00 P 1           1\n'.format(10.0e0*pbcbox[0]))
                    f.write(self.anum[1]+'\n')
                    for i,line in enumerate(self.lines[:-2]):
                        l = line[:30] + '{0:8.3f}{1:8.3f}{2:8.3f}'.format(10.0e0*pos[i][0], 10.0e0*pos[i][1], 10.0e0*pos[i][2])
                        l += line[56:]+'\n'
                        f.write(l)
                    f.write('TER\nENDMDL')

        if check_eneflag == True:
            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy()
            #print("thread name is "+str(threading.current_thread().name))
            #print('energy:', energy)
        else:
            energy = None

        log_f.close()
        if xtc_flag == True:
            xtc_reporter.close()
        elif hdf5_flag == True:
            hdf5_reporter.close() 
        # initialize simulaition.reportes
        simulation.reporters = []

        return simulation, energy


    def convert_gro2pdb(self, gro_):
        fname = gro_
        root, ext = os.path.splitext(gro_)
        outf = root + '.pdb'
        grostructure = pmd.load_file(gro_)
        grostructure.save(outf, format='PDB', overwrite=True)


    def getIntegratorInfo(self, integrator):
            dt = self.dt
            fric_const = integrator.getFriction()
            temperature = integrator.getTemperature()

            return dt, fric_const, temperature
