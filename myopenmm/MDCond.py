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
from numpy.random import *
import math

#  Requirement:
#  python 2.7, openmm, mdtraj, parmed
#

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


    def preparation(self, sysdir='SYS/', mddir='MD/', mode=None, index=''):
        
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
        if mode == None:
            sysgro = sysdir + self.sysfile
            systop = sysdir + self.forcefield
        elif mode == 'REMD':
            root, ext = self.sysfile.split('.')
            ext = '.' + ext
            sysgro = sysdir + root + index + ext
            systop = sysdir + self.forcefield

        # System output
        mddir = 'MD/'

        return sysgro, systop


    def create_top(self, systop, gro):
        if self.pbc:
            print('set periodic boundary condition...')
            top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            top = GromacsTopFile(systop)

        self.top = top
        return top


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
        gro = GromacsGroFile(sysgro)

        # Create top
        top = self.create_top(systop, gro)

        # Create system
        if self.heavyhydrogen:
            print('Hydrogen Mass Repartitioning...')
            self.hmass = 4*amu
        else:
            self.hmass = None

        if self.nonbonded_method == 'PME':
            print('set PME...')
            if self.nonbonded_cutoffflag:
                nonbonded_cutoff = self.nonbonded_cutoff   
                system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=PME, 
                                          nonbondedCutoff=self.nonbonded_cutoff,
                                          constraints=self.constraints)
            else:
                system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=PME,
                                          constraints=self.constraints)        
        elif self.nonbonded_method == 'Cutoff':
            print('set cutoff...')
            system = top.createSystem(hydrogenMass=self.hmass,nonbondedMethod=Cutoff, 
                                      nonbondedCutoff=self.nonbonded_cutoff,
                                      constraints=self.constraints)

        self.system = system

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
        simulation.minimizeEnergy(maxIterations=5000)

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

    def mdrun(self, simulation, ensname, index, mddir='MD/', sysdir='SYS/', 
              assert_system=True, check_eneflag=False, tstate=None):
        # ensname simulation
#        print('self.steps:', self.steps, 'mdresctep:', self.recstep)
        print(ensname+' Simulation...')
        mdname = ensname + index
        mdlog = mddir + mdname + '.log'
        log_f = open(mdlog, mode='a+')
        log_reporter = StateDataReporter(log_f, self.recstep, time=True,
                                                             totalEnergy=True, temperature=True, density=True,
                                                             progress=True, remainingTime=True, speed=True, totalSteps=self.steps, separator='\t')
        simulation.reporters.append(log_reporter)

        xtc_flag = False
        print('tstate:', tstate)
        if self.recflag and tstate == None:
            xtc_flag = True
            print('\nSaving...')
            mdxtc = mddir + mdname + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, self.recstep)
            simulation.reporters.append(xtc_reporter)
        elif self.recflag and tstate == 'SystemID:01':
            xtc_flag = True
            print('\nREMD Saving...')
            mdxtc = mddir + mdname + '.xtc'
            xtc_reporter = mymm.XTCReporter(mdxtc, self.recstep)
            simulation.reporters.append(xtc_reporter)
       
        if assert_system == True:
            print('Assert system...')
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
        del pdbstructure
        os.remove(mdpdb)

        if check_eneflag == True:
            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy()
        else:
            energy = None

        log_f.close()
        if xtc_flag == True:
            xtc_reporter.close()
        
        # initialize simulaition.reportes
        simulation.reporters = []
        return simulation, energy


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


class REMDConductor(MDConductor, object):
    def __init__(self, T_list):
        super(REMDConductor, self).__init__()
        self.Ts = list(map(float, T_list))
        self.n_replica = len(self.Ts)
        # Statistics
        self.attempts = [0 for i in range(self.n_replica)]
        self.successes = [0 for i in range(self.n_replica)]
        self.probability = [0.0 for i in range(self.n_replica)]

    def spread_replicas(self, simulation):
        print(self.n_replica, self.Ts)
        simulations = [ '' for i in range(self.n_replica)]

        dt, fric_const, temperature = self.getIntegratorInfo(simulation.integrator)
        print(dt, fric_const, temperature)
        system = simulation.system
        top = simulation.topology
        pos = simulation.context.getState(getPositions=True).getPositions()

        for i, T_ in enumerate(self.Ts):
            integrator_ = LangevinIntegrator(T_, fric_const, dt)
            properties = {'Precision': self.precision}
            simulations[i] = Simulation(top, system, integrator_, self.pltform,  properties)
            simulations[i].context.setPositions(pos)

        return simulations

       
    def make_arglist(self, simulations, ensname, index, mddir='MD/', sysdir='SYS/'): 
        arglist = [ []*5 for line in range(self.n_replica)]
        for i, T_ in enumerate(self.Ts):
            j = '{0:02d}'.format(i+1)
            indexj = index + '_' + j
            arglist[i] = [simulations[i], ensname, indexj,  mddir, sysdir]

        return arglist


    def initialize_replicas(self, simulations):
        tstates = [[' ', ' '] for i in range(self.n_replica)]
        for i, sim in enumerate(simulations):
            tstates[i][0] = 'SystemID:{0:02d}'.format(i+1)
            tstates[i][1] = sim
        return tstates 


    def remdrun(self, tstates, ensname, index, niter, mddir='MD', sysdir='MD'):
        simulations = [row[1] for row in tstates]
        arglist = self.make_arglist(simulations, ensname, index, mddir, sysdir)

        enes = []
        for i, args in enumerate(arglist):
            j = '{0:02d}'.format(i+1)
            k = '{0:04d}'.format(niter)
            indexj = index + '_' + j + '_' + k 
            sim, energy = self.mdrun(args[0], args[1], indexj, args[3], args[4], assert_system=False, check_eneflag=True, tstate=tstates[i][0])
            tstates[i][1] = sim
            enes.append(energy)

        return tstates, enes


    def del_intermediate(self, ensname, index, niter, niter_tot, mddir='MD/', sysdir= 'SYS/'):
        for j in range(self.n_replica):
            groname = sysdir + ensname + index + '_{0:02d}'.format(j+1) + '_{0:04d}'.format(niter) + '.gro'
            if not niter == niter_tot:
                os.remove(groname)


    def calc_prob(self, E_list, niter):
        k = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB
        b_list = [0 for i in range(self.n_replica)]
        p_list = [1 for i in range(self.n_replica)]
       
        if niter % 2 == 1:
            eo = 'odd'
        else:
            eo = 'even'

        if eo == 'even':
            elist = E_list[0::2]
        elif eo == 'odd':
            elist = E_list[1::2]
            elist = elist[:-1]

        for i,e in enumerate(E_list[:-1]):
            if eo == 'even' and i % 2 == 1:
                continue
            elif eo == 'odd' and i % 2 == 0:
                continue
            tm = self.Ts[i+1]*kelvin
            tn = self.Ts[i]*kelvin 
            Dbeta = (1/(k*tm) - 1/(k*tn))
            expbD = math.exp(Dbeta * (E_list[i+1] - E_list[i]))
            b_list[i] = expbD
            p = rand()
            p_list[i] = p
            
        return b_list, p_list

    def exchange(self, tstates, b_list, p_list, niter, index):
        if niter % 2 == 1:
            for i in range(1, self.n_replica, 2):
                self.attempts[i] += 1
        else:
            for i in range(0, self.n_replica, 2):
                self.attempts[i] += 1
        print('attempts:', self.attempts)

        simulations = [row[1] for row in tstates]
        exchange_flags = ['' for i in range(self.n_replica)]
        for i, b in enumerate(b_list):
            if b >= 1 or b >= p_list[i]:
                exchange_flags[i] = True
                self.successes[i] += 1
            else:
                exchange_flags[i] = False
        print('accepts:', self.successes)
        print('exchange_flags:', exchange_flags)

        # exchange position i-j position index
        positions = [ sim.context.getState(getPositions=True).getPositions() for sim in simulations]
        for i, flag in enumerate(exchange_flags):
            if flag == True:
                dummy_p = positions[i]
                dummy_i = tstates[i][0] 
                positions[i] = positions[i+1]
                tstates[i][0] = tstates[i+1][0]
                positions[i+1] = dummy_p
                tstates[i+1][0] = dummy_i

        for i,tstate in enumerate(tstates):
            tstate[1].context.setPositions(positions[i])
 
        # write exchange information
        states0 = '# st1\t\t' + '\t\t'.join(['st{0:d}'.format(i+1) for i in range(1, self.n_replica)]) 
        states = ['' for i in range(2*self.n_replica)]
        for j,i in enumerate([str(i+1) for i in range(self.n_replica)]):
            k = j + 1
            states[2*k-2] = i
            if j == self.n_replica + 1:
                break
            if exchange_flags[j]:
                states[2*k-1] = 'x'
                dummy = positions[j]
                positions[j] = positions[j+1]
                positions[j+1] = dummy

            elif exchange_flags[j] == False and j < self.n_replica - 1:
                states[2*k-1] = ' '

        states = '\t'.join(states) + '\n'
        exchange_fname = 'MD/exchange_'+index+'.log'
        if niter == 1:
            with open(exchange_fname, 'wt') as f:
                f.write(states0 + '\n')
        with open(exchange_fname, 'a+') as f:
            f.write(states)
 
        return tstates       

    def statistics(self):
        for i,s in enumerate(self.successes):
            self.probability[i] = '{0:4.3f}'.format(float(s)/float(self.attempts[i]))
        p = self.probability
        print('''##############################
       REMD  Statistics 
##############################
''')
        str_ = 'states:\t' + '\t'.join(['{0:02d}-{1:02d}'.format(i+1, i+2) for i in range(self.n_replica-1)])
        print(str_)
        str_p = 'prob:\t' + '\t'.join(self.probability[:-1])
        print(str_p)
        str_s = 'accept:\t' + '\t'.join(['{0:5d}'.format(s) for s in self.successes[:-1]])
        print(str_s)
        str_a = 'attmps:\t' + '\t'.join(['{0:5d}'.format(s) for s in  self.attempts[:-1]])
        print(str_a)


    def conduct(self, inpf, index, niter, mode='REMD', equilibriation=True):
        if mode == 'REMD':
            if equilibriation == True:
                self.loadFile(inpf)
                nvtgro, systop = self.preparation(sysdir='SYS/', mddir='MD/')
                simulation = self.setup(nvtgro, systop)
                self.mdrun(simulation, 'npt', index)
                simulations = self.spread_replicas(simulation)
            else:
                simulations = []
                self.loadFile(inpf)  
                for i in range(self.n_replica):
                   indexj = '{0:02d}'.format(i+1)
                   nvtgro , systop = self.preparation(sysdir='SYS/', mddir='MD/', mode=mode, index=indexj)
                   simulation = self.setup(nvtgro, systop)
                   simulations.append(simulation)
            tstates = self.initialize_replicas(simulations)

            for iter_ in range(1, niter+1):
                print('iter:', iter_)
                tstates, energys = self.remdrun(tstates, 'npt', index, iter_, mddir='MD/', sysdir='SYS/')
                b_list, p_list = self.calc_prob(energys, iter_)
                tstates = self.exchange(tstates, b_list, p_list, iter_, index)
                self.del_intermediate('npt', index, iter_, niter, mddir='MD/', sysdir='SYS/')

            # Check Statistics
            self.statistics()
