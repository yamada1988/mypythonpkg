import myopenmm as mymm
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import sys 

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
                        'platform' : None, 
                        'precision' : 'double'}

    def loadFile(self, file):
        try:
            f = open(file, 'rt')
        except IOError:
            sys.exit('Input file {0} not founded.'.format(file))

        with open(file, 'rt') as f:
            InpDict = {line.split('=')[0].strip() : line.split('=')[1].strip() \
                           for line in f if "#" not in line and "=" in line}

        InpDict.update(self.InpDict)
 
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
        self.nonbondedswitch = float(InpDict['nonbonded_switch'])
        self.reporterformat = InpDict['reporterformat']
        self.constraints = InpDict['constraints']
        self.nonbonded_method = InpDict['nonbonded_method']
        self.nonbonded_cutoffflag =InpDict['nonbonded_cutoffflag']
        self.nonbonded_cutoff = float(InpDict['nonbonded_cutoff'])
        self.heavyhydrogen = InpDict['heavyhydrogen']
        self.removeCMMotion = InpDict['removeCMMotion'] 
        self.platform = InpDict['platform']
        self.precision = InpDict['precision']

    def conduct(self, sysdir='SYS/', 
                emname='em', nvtname='nvt', nptname='npt', mdname='md', 
                mddir='MD/'):

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

        if self.mdstep:
            mdstep  = self.mdstep
        else:
            sys.exit('mdstep must be specified in input file.')

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
            top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            top = GromacsTopFile(systop)
       
        if self.heavyhydrogen:
            hmass = 4*amu
        else:
            hmass = 1*amu

        if self.nonbonded_method == 'PME':
            if self.nonbonded_cutoffflag:
                system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=PME, 
                                          nonbondedCutoff=self.nonbonded_cutoff*nanometer,
                                          constraints=self.constraints)
            else:
                system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=PME,
                                          constraints=self.constraints)        
        elif self.nonbonded_method == 'Cutoff':
            system = top.createSystem(hydrogenMass=hmass,nonbondedMethod=Cutoff, 
                                      nonbondedCutoff=self.nonbonded_cutoff*nanometer,
                                      constraints=self.constraints)


        if self.integrator == 'Langevin':
            integrator = LangevinIntegrator(temperature, fric_const, dt) 
        elif self.integrator == 'Brownian':
            integrator = BrownianIntegrator(temperature, fric_const, dt)
        elif self.integrator == 'Velret':
            integrator = VelretIntegrator(dt)
        else:
            sys.exit('Invalid Integrator type. Check your input file.')

        # EM simulation
        if self.emflag:
            simulation = Simulation(top.topology, system, integrator)
            simulation.context.setPositions(gro.positions)

            print('Minimizing...')
            empdb = mddir + 'em.pdb'
            simulation.minimizeEnergy()

            print('Saving...')
            positions = simulation.context.getState(getPositions=True).getPositions()
            PDBFile.writeFile(simulation.topology, positions, open(empdb, 'w'))

        # NVT simulation
        if self.nvtflag:
            print ('NVT Simulation...')
            if not self.emflag:
                simulation = Simulation(top.topology, system, integrator)
                simulation.context.setPositions(gro.positions)
            else:
                simulation.context.setVelocitiesToTemperature(temperature)
             
            if self.nvtrecflag:
                nvtpdb = mddir + nvtname + '.pdb'
                nvtlog = mddir + nvtname + '.log'
                simulation.reporters.append(StateDataReporter(nvtlog, nvteqrecstep, step=True,
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
            simulation = Simulation(top.topology, system, integrator)


            if not self.emflag and not self.nvtflag:
                positions = gro.positions

            simulation.context.setPositions(positions)

            if self.nptrecflag:
                nptpdb = mddir + nptname + '.pdb'
                nptlog = mddir + nptname + '.log'
                simulation.reporters.append(StateDataReporter(nptlog, npteqrecstep, step=True,
                                                              totalEnergy=True, temperature=True, density=True, 
                                                              progress=True, remainingTime=True, speed=True, 
                                                              totalSteps=nptstep, separator='\t'))
            nptchk = mddir + nptname + '.chk'
            simulation.reporters.append(CheckpointReporter(nptchk, npteqrecstep))

            if self.nptrecflag:
                print('\nSaving...')
                nptxtc = mddir + nptname + '.xtc'
                xtc_reporter = mymm.XTCReporter(nptxtc, npteqrecstep)
                simulation.reporters.append(xtc_reporter)

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

            simulation = Simulation(top.topology, system, integrator)
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
 
