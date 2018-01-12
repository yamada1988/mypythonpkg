import myopenmm as mymm
import mdtraj as md
from distutils.util import strtobool
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import sys 

class MDConductor_legacy:
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
                        'pme-alpha' : 0.0, 
                        'pme-nx' : 64,
                        'pme-ny' : 64,
                        'pme-nz' : 64,
                        'pme-etol' : 1.0e-5,     
                        'nvtstep' : 5000, 
                        'nptstep' :50000, 
                        'mdstep' : 50000, 
                        'nvteqrecstep' : 500, 
                        'npteqrecstep' : 500,
                        'mdrecstep' : 500, 
                        'emflag' : True, 
                        'nvtflag' : True, 
                        'nptflag' : True, 
                        'nvtrecflag' : True, 
                        'nptrecflag' : True, 
                        'mdrecflag' :  True, 
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
                        'precision' : 'double'}
        
        print('Initialize...')               
    def loadFile(self, file):
        try:
            f = open(file, 'rt')
        except IOError:
            sys.exit('Input file {0} not founded.'.format(file))

        with open(file, 'rt') as f:
            InpDict = {line.split('=')[0].strip() : line.split('=')[1].strip() \
                           for line in f if "#" not in line and "=" in line}

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
        self.pme_alpha = float(InpDict['pme-alpha'])
        self.pme_nx = int(InpDict['pme-nx'])
        self.pme_ny = int(InpDict['pme-ny'])
        self.pme_nz = int(InpDict['pme-nz'])
        self.nvteqrecstep = int(InpDict['nvteqrecstep'])
        self.npteqrecstep = int(InpDict['npteqrecstep'])
        self.mdrecstep = int(InpDict['mdrecstep'])
        self.emflag = strtobool(InpDict['emflag'])
        self.nvtflag = strtobool(InpDict['nvtflag'])
        self.nptflag = strtobool(InpDict['nptflag'])
        self.nvtrecflag = strtobool(InpDict['nvtrecflag'])
        self.nptrecflag = strtobool(InpDict['nptrecflag'])
        self.mdrecflag = strtobool(InpDict['mdrecflag'])
        self.simlogflag = strtobool(InpDict['simlogflag'])
        self.simlogstep = int(InpDict['simlogstep'])
        self.verbos = InpDict['verbose']
        self.forcerecflag = strtobool(InpDict['forcerecflag'])
        self.pbc = InpDict['pbc']
        self.remdflag = strtobool(InpDict['remdflag'])
        self.annealingflag = strtobool(InpDict['annealingflag'])
        self.switchflag = strtobool(InpDict['nonbonded_switchflag'])
        self.nonbonded_switch = float(InpDict['nonbonded_switch'])
        self.reporterformat = InpDict['reporterformat']
        self.constraints = InpDict['constraints']
        self.nonbonded_method = InpDict['nonbonded_method']
        self.nonbonded_cutoffflag = InpDict['nonbonded_cutoffflag']
        self.nonbonded_cutoff = float(InpDict['nonbonded_cutoff'])
        self.heavyhydrogen = strtobool(InpDict['heavyhydrogen'])
        self.removeCMMotion = InpDict['removeCMMotion'] 
        self.ghost_particle = strtobool(InpDict['ghost_particle'])
        self.path_ndxfile = InpDict['path_ndxfile']
        self.platform = InpDict['platform']
        self.precision = InpDict['precision']

        with open('openmm.out', 'wt') as f:
            print('writing openmm.out...')
            for k, v in InpDict.items():
                k = k.ljust(20)
                f.write('{0}:{1}\n'.format(k, v))


    def conduct(self, sysdir='SYS/', 
                emname='em', nvtname='nvt', nptname='npt', mdname='md', 
                mddir='MD/'):

        # Platform input
        pltfmname = self.platform
        precision = self.precision
        platform = Platform.getPlatformByName(pltfmname)
        if self.platform == 'CUDA':
            properties = {'CudaPrecision': precision}#, 'CudaDeviceIndex': '0,1,2,3'}#, 'CudaUseBlockingSync':False, 'CudaDeviceIndex':'0,1,2,3'}
        elif self.platform == 'CPU':
            properties = {}
        elif self.platform == 'OpenCL':
            properties = {'OpenCLPrecision': precision}#, 'OpenCLDeviceIndex': '0,1,2,3'}


        # System input
        sysgro = sysdir + self.sysfile
        systop = sysdir + self.forcefield

        # System output
        mddir = 'MD/'

        # Simulation setting
        if self.temperature:
            temperature = self.temperature * kelvin

        if self.friction_const:
            fric_const = self.friction_const / picosecond

        dt = self.dt * picosecond

        if self.pressure:
            pressure = self.pressure * bar

        if self.nvtflag:
            nvtstep = self.nvtstep  
        else:
            nvtstep = 0 

        if self.nptflag:
            nptstep = self.nptstep
        else: 
            nptstep = 0
        
        mdstep  = self.mdstep

        totstep = nvtstep + nptstep + mdstep 

        if self.nvtrecflag:
            nvteqrecstep = self.nvteqrecstep
        if self.nptrecflag:
            npteqrecstep = self.npteqrecstep
        if self.mdrecflag:
            mdrecstep =  self.mdrecstep

        if self.simlogflag:
            simlogstep = self.simlogstep


        # Simulation Setting
        gro = GromacsGroFile(sysgro)
        if self.pbc:
            print('set periodic boundary condition...')
            top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            top = GromacsTopFile(systop)

       
        if self.heavyhydrogen:
            print('Hydrogen Mass Repartitioning...')
            hmass = 4*amu
        else:
            hmass = 1*amu


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


        if self.switchflag:
            print('set switching functio...')
            nonbonded_switch = self.nonbonded_switch
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            forces['NonbondedForce'].setUseSwitchingFunction(True)
            forces['NonbondedForce'].setSwitchingDistance(nonbonded_switch)
 

        if self.integrator == 'Langevin':
            print('set Langevin Integrator...')
            integrator = LangevinIntegrator(temperature, fric_const, dt) 
        elif self.integrator == 'Brownian':
            print('set Brownian Integrator...')
            integrator = BrownianIntegrator(temperature, fric_const, dt)
        elif self.integrator == 'Velret':
            print('set Verlet Integrator...')
            integrator = VelretIntegrator(dt)
        else:
            sys.exit('Invalid Integrator type. Check your input file.')

        # Ghost-particle        
        if self.ghost_particle:
            print('set ghost particle calculation...')
            inpf = self.path_ndxfile        
            with open(inpf, 'rt') as gf:
                total_lines =[line.strip() for line in gf]
                try:
                    ghost_index = total_lines.index('[ Ghost ]')
                    solvent_index = total_lines.index('[ Solvent ]')
                except:
                    sys.exit('If ghost_particle flag is True, ghost_index file must be specified.\n')
            ghost_particles = [int(line)-1 for line in ' '.join(total_lines[ghost_index+1:solvent_index]).split() ]
            solvent_particles = [int(line)-1 for line in ' '.join(total_lines[solvent_index+1:]).split() ]
#            print(ghost_particles)
#            print(solvent_particles)

            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            nonbonded_force = forces['NonbondedForce']

#            print('Before:')
#            print(nonbonded_force.getNumExceptions())
#            for i in range(nonbonded_force.getNumExceptions()):
#                print(nonbonded_force.getExceptionParameters(i))

            for sp in solvent_particles:
                for gp in ghost_particles:
                    nonbonded_force.addException(sp, gp, 0.0e0*0.0e0, sigma=1.0e0, epsilon=0.0e0)

#            print('After:')
#            print(system.getForce(0).getNumExceptions())

#            for i in range(nonbonded_force.getNumExceptions()):
#                print(nonbonded_force.getExceptionParameters(i))
                       

        # EM simulation
        if self.emflag:
            simulation = Simulation(top.topology, system, integrator, platform, properties)
            simulation.context.setPositions(gro.positions)

            print('Minimizing...')
            empdb = mddir + 'em.pdb'
            simulation.minimizeEnergy(maxIterations=20)

            print('Saving...')
            positions = simulation.context.getState(getPositions=True).getPositions()
            PDBFile.writeFile(simulation.topology, positions, open(empdb, 'w'))

        # NVT simulation
        if self.nvtflag:
            print ('NVT Simulation...')
            if not self.emflag:
                simulation = Simulation(top.topology, system, integrator, platform)
                simulation.context.setPositions(gro.positions)
            else:
                simulation.context.setVelocitiesToTemperature(temperature)
             
            if self.nvtrecflag:
                nvtpdb = mddir + nvtname + '.pdb'
                nvtlog = mddir + nvtname + '.log'
                simulation.reporters.append(StateDataReporter(nvtlog, nvteqrecstep, time=True,
                                                              totalEnergy=True, temperature=True, density=True, 
                                                              progress=True, remainingTime=True, speed=True, 
                                                              totalSteps=nvtstep, separator='\t'))

            simulation.step(nvtstep)

            if self.nvtrecflag:
                print('Saving...')
                positions = simulation.context.getState(getPositions=True).getPositions()
                PDBFile.writeFile(simulation.topology, positions, open(nvtpdb, 'w'))
 
        # NPT simulation
        if self.nptflag:
            print('NPT Simulation...')
            system.addForce(MonteCarloBarostat(pressure, temperature))
            if self.integrator == 'Langevin':
                integrator = LangevinIntegrator(temperature, fric_const, dt) 
            elif self.integrator == 'Brownian':
                integrator = BrownianIntegrator(temperature, fric_const, dt)
            elif self.integrator == 'Velret':
                integrator = VelretIntegrator(dt)
            else:
                sys.exit('Invalid Integrator type. Check your input file.')
            simulation = Simulation(top.topology, system, integrator, platform, properties)


            if not self.emflag and not self.nvtflag:
                positions = gro.positions

            simulation.context.setPositions(positions)

            if self.nptrecflag:
                nptlog = mddir + nptname + '.log'
                simulation.reporters.append(StateDataReporter(nptlog, npteqrecstep, time=True,
                                                              totalEnergy=True, temperature=True, density=True, 
                                                              progress=True, remainingTime=True, speed=True, 
                                                              totalSteps=nptstep, separator='\t'))
            nptchk = mddir + nptname + '.chk'
            simulation.reporters.append(CheckpointReporter(nptchk, npteqrecstep))

            if self.nptrecflag:
                print('\nSaving...')
                if self.reporterformat == 'XTC':
                    print('Save as XTC file')
                    nptxtc = mddir + nptname + '.xtc'
                    xtc_reporter = mymm.XTCReporter(nptxtc, npteqrecstep)
                    simulation.reporters.append(xtc_reporter)
                elif self.reporterformat == 'PDB':
                    print('Save as PDB file')
                    nptpdb = mddir + nptname + '.pdb'
                    pdb_reporter = PDBReporter(nptpdb, npteqrecstep) 
                    simulation.reporters.append(pdb_reporter)

            simulation.step(nptstep)

        # MD simulation
        print('MD Simulation...')
        if not self.emflag and not self.nvtflag and not self.nptflag:
            if mdensemble == 'npt':
                system.addForce(MonteCarloBarostat(pressure, temperature))
                if self.integrator == 'Langevin':
                    integrator = LangevinIntegrator(temperature, fric_const, dt) 
                elif self.integrator == 'Brownian':
                    integrator = BrownianIntegrator(temperature, fric_const, dt)
                elif self.integrator == 'Velret':
                    integrator = VelretIntegrator(dt)
                else:
                    sys.exit('Invalid Integrator type. Check your input file.')

            simulation = Simulation(top.topology, system, integrator, platform)
            simulation.context.setPositions(gro.positions)
        
        if self.simlogflag:
            mdlog = mddir + mdname + '.log'
            simulation.reporters.append(StateDataReporter(mdlog, mdrecstep, time=True,
                                                          totalEnergy=True, temperature=True, density=True, 
                                                          progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))

        mdchk = mddir + mdname + '.chk'
        simulation.reporters.append(CheckpointReporter(mdchk, mdrecstep))

        if self.mdrecflag:
            print('\nSaving...')
            mdxtc = mddir + mdname + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, mdrecstep)
            simulation.reporters.append(xtc_reporter)

        simulation.step(mdstep)
        print('Done!\n')  
 
