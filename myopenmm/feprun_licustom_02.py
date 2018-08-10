import myopenmm as mymm
import math
import mdtraj as md
from distutils.util import strtobool
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from sys import stdout
import os
import sys 

class FEPConductor_ljcustom:
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
                        'fepflag': False,
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
        self.setPME = strtobool(InpDict['setPME'])
        self.pme_alpha = float(InpDict['pme-alpha'])
        self.pme_nx = int(InpDict['pme-nx'])
        self.pme_ny = int(InpDict['pme-ny'])
        self.pme_nz = int(InpDict['pme-nz'])
        self.pme_ftol = float(InpDict['pme-ftol'])
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
        self.pbc = strtobool(InpDict['pbc'])
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
        self.fepflag = strtobool(InpDict['fepflag'])
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
                mddir='MD/', statenum=0):
        
        print('conduct start...')
        # Platform input
        pltfmname = self.platform
        precision = self.precision
        platform = Platform.getPlatformByName(pltfmname)
        if self.platform == 'CUDA':
            properties = {'CudaPrecision': precision, 'CudaDeviceIndex': '0,1'}#, 'CudaUseBlockingSync':False }
        elif self.platform == 'CPU':
            properties = {}
        elif self.platform == 'OpenCL':
            properties = {'OpenCLPrecision': precision}#, 'OpenCLDeviceIndex': '0,1,2,3'}


        # System input
        sysgro = sysdir + self.sysfile
        systop = sysdir + self.forcefield

        # System outputdir
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
        if self.fileformat == 'GRO':
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
            hmass = None


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


        if self.integrator == 'Langevin':
            print('set Langevin Integrator...')
            integrator0 = LangevinIntegrator(temperature, fric_const, 0.00200) 
        elif self.integrator == 'Brownian':
            print('set Brownian Integrator...')
            integrator0 = BrownianIntegrator(temperature, fric_const, 0.00200)
        elif self.integrator == 'Verlet':
            print('set Verlet Integrator...')
            integrator0 = VerletIntegrator(0.00200)
        else:
            sys.exit('Invalid Integrator type. Check your input file.')


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

            if int(ghost_start) == 0 and int(ghost_start) == int(ghost_end):
                ghost_particles = range(0,0)
            else:
                ghost_particles = range(int(ghost_start)-1, int(ghost_end))
            solvent_particles = range(int(solvent_start)-1, int(solvent_end))

#            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
#            nonbonded_force = forces['NonbondedForce']

#            solute_particles = set.union(set(core_particles), set(ghost_particles))
#            slt_param = [[]]*len(solute_particles) 
#            for index in solute_particles:
#                slt_param[index] = nonbonded_force.getParticleParameters(index)
#                if index in ghost_particles:
#                    nonbonded_force.setParticleParameters(index, charge=0.0e0, sigma=slt_param[index][1], epsilon=0.0e0)
#
#            for i in ghost_particles:
#                for j in solute_particles:
#                    if i == j:

#                        continue
#                    try:
#                        nonbonded_force.addException(i, j, chargeProd=slt_param[i][0]*slt_param[j][0], 
#                                                           sigma=0.50*(slt_param[i][1]+slt_param[j][1]),
#                                                           epsilon=sqrt(slt_param[i][2]*slt_param[j][2]))
#                    except:
#                        pass

        # Retrieve the NonbondedForce
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        nbforce = forces['NonbondedForce']

        # Add a CustomNonbondedForce to handle only alchemically-modified interactions
        alchemical_particles = set([0,1,2,3,4,5,6,7,8])
        print('alchemical particles:', alchemical_particles)
        chemical_particles = set(range(system.getNumParticles())) - alchemical_particles
        energy_function = '0.0e0*charge1*charge2*0.0e0 + lambda_vdw * 4 * epsilon * x * (x-1.0); x = ((sigma^(6.0e0))/reff_sterics);'
        energy_function += 'reff_sterics = sigma^(6.0e0)*0.50e0*(1-lambda_vdw)+r^(6.0e0);'
        energy_function += 'sigma = 0.50e0*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
        custom_force = openmm.CustomNonbondedForce(energy_function)

        custom_force.addGlobalParameter('lambda_vdw', 1.0e0)
        custom_force.addPerParticleParameter('charge')
        custom_force.addPerParticleParameter('sigma')
        custom_force.addPerParticleParameter('epsilon')
        nb_params = []*system.getNumParticles()
        for index in range(system.getNumParticles()):
            nb_params.append(nbforce.getParticleParameters(index))
            custom_force.addParticle([nb_params[index][0]*elementary_charge, nb_params[index][1]*nanometer, nb_params[index][2]*kilojoule/mole])
        custom_force.addInteractionGroup(alchemical_particles, chemical_particles)

        # Retrieve NonbondedForce beyond 1-4 interactions
        print('Set interaction parameters beyond 1-4 interactions in alchemical particles...')
        for i in alchemical_particles:
            for j in alchemical_particles:
                if i == j:
                    continue
                try:
                    nbforce.addException(i, j, chargeProd=nb_params[i][0]*nb_params[j][0],
                                         sigma=0.50*(nb_params[i][1]+nb_params[j][1]),
                                         epsilon=sqrt(nb_params[i][2]*nb_params[j][2]))
                    print('pair:', i, j, nb_params[i], nb_params[j])
                except:
                    pass

        print('Set NonbondedMethod for Modified-LJ interaction')
        print('before:{0:d}'.format(custom_force.getNonbondedMethod()))
        custom_force.setNonbondedMethod(2)
        print('after:{0:d}'.format(custom_force.getNonbondedMethod()))
        custom_force.setCutoffDistance(1.20e0 * nanometer)
        custom_force.setUseSwitchingFunction(True)
        custom_force.setSwitchingDistance(1.0e0 * nanometer)#nonbonded_switch)

        print('Check exception parameters for nbforce...')
        for index in range(nbforce.getNumExceptions()):
            exception_info = nbforce.getExceptionParameters(index)
            print(exception_info)

        # Set Exculusion Info for Modified LJ and Coulomb interactions from NonbondedForce
        print('Check exclusion info for NonbondedForce')
        for index in range(nbforce.getNumExceptions()):
            exception_info = nbforce.getExceptionParameters(index)
            index_v = custom_force.addExclusion(exception_info[0], exception_info[1])
            custom_force.setExclusionParticles(index_v, exception_info[0], exception_info[1])
            print(index_v, exception_info[0], exception_info[1])

        system.addForce(custom_force)


        # Collect data 
        kmax = 1
        equiliter = 100 # number of steps per state
        nsteps = 2000 # number of steps per sample
        niterations = 301 # number of samples to collect per alchemical state
        lambdas_q = np.linspace(1.0e0, 0.0e0, 11) 
        lambdas_v = np.linspace(1.0e0, 0.0e0, 11) 
        nlamq = len(lambdas_q)
        nlamv = len(lambdas_v)
        print('lambdas_q:', lambdas_q)
        print('lambdas_v:', lambdas_v)
        nstates = len(lambdas_q) + len(lambdas_v)
        u_kln = np.zeros([nstates,nstates,niterations], np.float64)
        kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator0.getTemperature()
        snum = int(statenum)

        # Add a barostat
        barostat = openmm.MonteCarloBarostat(pressure, temperature)
        system.addForce(barostat)

        # Create a contextq
        simulation0 = Simulation(top.topology, system, integrator0, platform, properties)
        simulation0.context.setPositions(gro.positions)
        simulation0.context.setVelocitiesToTemperature(temperature)

        print('Minimizing energy...')
        simulation0.minimizeEnergy(tolerance = 1.0, maxIterations=5000)

        # Conduct simulation in forward direction
        for k in range(snum, snum+1):
            stored_fname = 'MD/u_{0:03d}.xvg'.format(k)
            with open(stored_fname, 'wt') as f:
                f.write('#state1\tstate2\titer\tu_kln (kT)\n')

            for iteration in range(niterations):
                if k < nlamq - 1:
                    lambda_q = lambdas_q[k] 
                    lambda_v = 1.0e0
                elif k >= nlamq - 1:
                    lambda_q = 0.0e0
                    lambda_v = lambdas_v[k-nlamq+1]

                print('Set alchemical state in NPT simulation...')
                # Retrieve the NonbondedForce
                print('Set lamda-modified parameters for LJ interaction in NPT simulation...')
                forces = { force.__class__.__name__ : force for force in system.getForces() }
                nbforce = forces['NonbondedForce']
                for index in range(system.getNumParticles()):
                    if index in alchemical_particles:
                        nbforce.setParticleParameters(index, lambda_q*nb_params[index][0]*elementary_charge, nb_params[index][1]*nanometer, 0.0e0*kilojoule/mole)
                        print(nbforce.getParticleParameters(index))
                    else:
                         nbforce.setParticleParameters(index, nb_params[index][0]*elementary_charge, nb_params[index][1]*nanometer, nb_params[index][2]*kilojoule/mole)    

                # Update Parameter value in Context
                context_0 = simulation0.context
                context_0.setParameter('lambda_vdw', lambda_v)
                nbforce.updateParametersInContext(context_0)

                print('state {0:02d} iteration {1:02d} / {2:02d}\tlambda_q:{3:4.2f}\tlambda_v:{4:4.2f}'.format(k, iteration, niterations, lambda_q, lambda_v))

                mdlog = mddir + 'md{0:03d}'.format(k,) + '.log'

                simulation0.reporters.append(StateDataReporter(mdlog, 250, time=True,
                                                               potentialEnergy=True, kineticEnergy=True, temperature=True, density=True,
                                                               progress=True, remainingTime=True, speed=True,
                                                               totalSteps=nsteps*niterations, separator='\t'))

#                if k == 0 or k == nstates-2:
                if k == 0:
                    print('\nSaving...')
                    mdxtc = mddir + mdname + '_{0:02d}'.format(k) + '_{0:03d}'.format(iteration) + '.xtc'
                    xtc_reporter = mymm.XTCReporter(mdxtc, 50)
                    simulation0.reporters.append(xtc_reporter)

                # Run some dynamics
                simulation0.step(nsteps)

                # Check nbforce parameter
                nb_params_test = []*system.getNumParticles()
                for index in range(system.getNumParticles()):
                    nb_params_test.append(nbforce.getParticleParameters(index))
                    if index in alchemical_particles:
                        print(k, index, nb_params_test[index][0])

                # Compute energies at all alchemical states
                if iteration < equiliter:
                    pass

                elif iteration >= equiliter:
                    for l in range(nstates-1):
                        if l < nlamq - 1:
                            lambda_ql = lambdas_q[l]
                            lambda_vl = 1.0e0
                        elif l >= nlamq - 1 :
                            lambda_ql = 0.0e0 
                            lambda_vl = lambdas_v[l-nlamq+1]
                        context_l = simulation0.context
                        context_l.setParameter('lambda_vdw', lambda_vl)
                     
                        print('Change Parameters in nbforce for each state l in NPT simulation...')
                        print('state:{0:02d}'.format(l))
                        for index in range(system.getNumParticles()):
                            if index in alchemical_particles:
                                nbforce.setParticleParameters(index, lambda_ql*nb_params[index][0]*elementary_charge, nb_params[index][1], 0.0e0*kilojoule/mole)
                                print(nbforce.getParticleParameters(index))
                            else:
                                nbforce.setParticleParameters(index, nb_params[index][0]*elementary_charge, nb_params[index][1]*nanometer, nb_params[index][2]*kilojoule/mole)    


                        # Update Parameter value in Context
                        nbforce.updateParametersInContext(context_l)

                        u_kln[k,l,iteration] = context_l.getState(getEnergy=True).getPotentialEnergy() / kT
                        with open(stored_fname, 'a+') as f:
                            f.write('{0: 2d}\t{1:2d}\t{2:4d}\t{3:12.7f}\n'.format(k, l, iteration-equiliter, u_kln[k,l,iteration]))

                        del context_l

                del context_0 
