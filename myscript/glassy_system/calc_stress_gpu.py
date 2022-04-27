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

# ============================================================================
# Input filename
# ============================================================================
args = sys.argv
inpf = args[1]

infh = h5py.File(inpf, 'r')

# number of total particles            
nparticles = len(infh['coordinates'][0])
print(""" ###########################
   WAHNSTROM BINARY SYSTEM 
   Num: {0:9d}
 ###########################""".format(nparticles))

# =============================================================================
# Specify simulation parameters
# =============================================================================

n1 = int(nparticles/2)
n2 = n1

# species1
mass1 = 39.9 * amu # mass
sigma1 =0.34 * nanometers # Lennard-Jones sigma
epsilon1 = 0.994 * kilojoule/mole # Lennard-Jones well-depth
charge1 = 0.0e0 * elementary_charge # argon model has no charge

# species2
mass2 = 79.896 * amu # mass
sigma2 = 0.408 * nanometers # Lennard-Jones sigma
epsilon2 = 0.994 * kilojoule/mole# Lennard-Jones well-depth
charge2 = 0.0e0 * elementary_charge # argon model has no charge
sigma = cp.asarray(np.concatenate([(sigma1/nanometers)*np.ones(n1), (sigma2/nanometers)*np.ones(n2)]), dtype=np.float32)
epsilon = epsilon1/ (kilojoule/mole)
print('ep:',epsilon)
masses = np.concatenate([(mass1/amu)*np.ones(n1), (mass2/amu)*np.ones(n2)])

# temperature
temp = float(args[1].split('_')[1].split('.')[0])

# Argon units
kB = BOLTZMANN_CONSTANT_kB
temp0 = epsilon1 / (kB*AVOGADRO_CONSTANT_NA)
NA0 = AVOGADRO_CONSTANT_NA / mole
m0 = mass1
l0 = sigma1
e0 = epsilon1
Tstar = float(temp)
print("Target temperature:", Tstar)
kBT = kB * Tstar * kelvin * AVOGADRO_CONSTANT_NA / (kilojoule/mole) 

temperature = Tstar * kelvin # temperature
collision_rate = 2.0 / picosecond # collision rate for Langevin thermostat
timestep = 2.0 * femtosecond # integrator timestep

cutoff = 1.0 # Compute cutoff
print("sigma1, sigma2, cutoff")
print(sigma1, sigma2, cutoff)
lmax = 30

fhstep  =  20
atmstep =   1
total_steps = len(infh['coordinates'])
t_ratio = fhstep/atmstep
units_ = {'length': 'nanometers',
          'mass': 'amu',
          'time': 'picosecond',
          'force': 'kilojoule/mole/nanometer',
          'energy': 'kilojoule/mole',
          'dt_fh/dt_atm': t_ratio}
dt_atom = atmstep*timestep  / femtosecond

L0 = infh['cell_lengths'][0,0]
rN = int(L0/cutoff)
n_mesh = rN**3
outf = 'MD/test_local0001_{0:03d}.h5'.format(int(temp))
data = cp.arange(n_mesh, dtype=np.int32)
with h5py.File(outf, 'w') as f:
    f.create_group('Units')
    for k,v in units_.items():
        f['Units'].create_dataset(k, data=v)
    f.create_dataset('dt_atom', data=dt_atom)
    f.create_dataset('localMass', shape=(total_steps/fhstep, rN, rN, rN))
    f.create_dataset('localStress_v', shape=(total_steps/fhstep, rN, rN, rN, 3, 3))
    f.create_dataset('localStress_f', shape=(total_steps/fhstep, rN, rN, rN, 3, 3))
    f.create_dataset('localMomentum', shape=(total_steps/fhstep, rN, rN, rN, 3))
    f.create_dataset('MeshGrid', shape=(3, rN+1))
    f.create_dataset('time', shape=(total_steps/fhstep,))

#for it in range(nsample_steps/fhstep):
for it in range(total_steps/fhstep):
    print('it:',it)
    time_ = (it+1)*fhstep 
    lm = np.zeros((rN,rN,rN))
    lg = np.zeros((rN,rN,rN, 3))
    local_mass = 0.0e0
    local_g = 0.0e0
    local_v = 0.0e0
    local_sigma0k = np.zeros((rN, rN, rN, 3, 3))
    local_sigmav = np.zeros((rN, rN, rN, 3, 3))
    for t in range(fhstep/atmstep):
        print('t:',t)
        t_ = fhstep/atmstep * it + t
        pos = infh['coordinates'][t_]
        vel = infh['velocities'][t_]
        box = infh['cell_lengths']
        L = box[0,0]

        rN = int(L/cutoff)
        # unset PBC 
        pos -= np.trunc(pos/L)*L
        pos += np.round((0.50e0*L-pos)/L)*L
        pos = cp.asarray(pos, dtype=np.float32)

        # grid parameters
        r_min = 0.0e0 
        r_max = L
        dr = (r_max-r_min)/ float(rN)
        dV = dr**3
        if t == 0:
            r_ = np.array([r_min+dr*ir for ir in range(rN)])
            rax = np.array([r_, r_, r_])
        # calculate rho, g, v, sigma_ab, tau_ab, S_ab, eta as local variables
        m = cp.asarray(masses, dtype=np.float32)
        lm = cp.zeros((rN,rN,rN), dtype=cp.float32)
        r0 = cp.asarray(pos[:,0]/dr, dtype=cp.int32)
        r1 = cp.asarray(pos[:,1]/dr, dtype=cp.int32)
        r2 = cp.asarray(pos[:,2]/dr, dtype=cp.int32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, raw float32 m',
            'raw float32 lm',
            '''
              int ind[] = {r0,r1,r2};
              atomicAdd(&lm[ind], m[i]);
            '''            
            )(r0,r1,r2,m,lm)
        lm = cp.asnumpy(lm)
        lm /= dV
        local_mass += lm
        del lm

        v = cp.asarray(vel, dtype=np.float32)
        lg = cp.zeros((rN,rN,rN,3), dtype=cp.float32)
        cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, raw float32 m, raw float32 v',
            'raw float32 lg',
            '''
            for (int j = 0; j < 3; j++) {
              int ind[] = {r0, r1, r2, j};
              int ind2[2] = {i,j};
              atomicAdd(&lg[ind], m[i]*v[ind2]);
            }
            '''
            )(r0,r1,r2, m, v,lg)
        lg = cp.asnumpy(lg)

        local_g += lg /dV
        del lg

        ls0k = cp.zeros((rN,rN,rN,3,3), dtype=cp.float32)
        cp.ElementwiseKernel(
                  'int32 r0, int32 r1, int32 r2, raw float32 m, raw float32 v',
                  'raw float32 ls0k',
                  '''
                  for (int xa=0; xa<3; xa++){
                    for (int xb=0; xb<3; xb++){
                      int ind[] = {r0,r1,r2, xa, xb};
                      atomicAdd(&ls0k[ind], m[i]*v[3*i+xa]*v[3*i+xb]);
                    }
                  }
                  '''
                 )(r0,r1,r2, m, v, ls0k)
        ls0k = cp.asnumpy(ls0k)
        local_sigma0k += ls0k/dV 
        del ls0k

        address_list = -1*cp.ones((n_mesh, lmax), dtype=np.int32)
        cp.ElementwiseKernel(
            in_params='int32 data, raw int32 r0, raw int32 r1, raw int32 r2, float32 cutoff, int32 rN, int32 nparticles, int32 lmax', 
            out_params='raw int32 address_list',
            operation=\
            '''
            int l = 0;
            for (int j = 0; j < nparticles; j++){
              if ( i == r0[j] + rN*r1[j] + rN*rN*r2[j] ) {
                int ind[2] = {i, l};
                address_list[ind] = j;
                l += 1;
                } else {
                  continue;
                }
            } 
            ''',
            name='update_list')(data, r0, r1, r2, cutoff, rN, nparticles, lmax, address_list)

        cp.cuda.Stream.null.synchronize()

        stress = cp.zeros((n_mesh, 3, 3), dtype=np.float32)
        cp.ElementwiseKernel(
              in_params='int32 data, raw int32 address_list, int32 lmax, raw float32 sigma, float32 cutoff, float32 L, float32 epsilon, raw float32 pos, int32 rN', 
              out_params='raw float32 stress',
              operation=\
              '''
              int n_mesh = _ind.size();
              int pair_index[27] = {0.0};
              int i_count = 0;
              for (int ix0 = -1; ix0<2; ix0++){
                for (int iy0 = -1; iy0<2; iy0++){
                  for (int iz0 = -1; iz0<2; iz0++){
                    pair_index[i_count] = (i + rN*rN*iz0 +rN*iy0 + ix0)%n_mesh;
                    i_count += 1;
                  }
                }
              }
              for (int p1=0; p1<lmax; p1++){
                int ip = address_list[lmax*i+p1];
                if ( ip < 0 )
                {
                break;
                } 
                else 
                {
                  for (int p2=0; p2<lmax; p2++){
                    for (int j=0; j<27; j++){
                      int jp = address_list[lmax*pair_index[j]+p2];
                      if ( jp<0 || jp <= ip )
                      {
                      break;
                      }
                      else
                      {
                      double tmp[3] = {0.0};
                      double diff = 0.0;
                      for (int m = 0; m<3; m++){
                        tmp[m] = pos[3*ip+m] - pos[3*jp+m];
                        tmp[m] += L*int(tmp[m]/(0.50*L));
                        diff += powf(tmp[m],2); 
                      }
                      double r = sqrt(diff);
                      if ( r > cutoff)
                      {
                      continue;
                      }
                      else
                      {
                        double r6 = powf(0.50*(sigma[ip]+sigma[jp])/r, 6);
                        double f[3] = {0.0};
                        for (int ifx=0; ifx<3; ifx++){
                          f[ifx] = 24.0*epsilon*r6*(2*r6-1.0)/r*tmp[ifx]/r;
                          if ( f[ifx] > 50.0){
                            printf("%d, %f,   %f   :",i,  r/(0.50*(sigma[ip]+sigma[jp])), f[ifx]);
                          }
                        }
                        for (int xa=0; xa<3; xa++){
                          for (int xb=0; xb<3; xb++){
                            int ind_mesh1[3] = {i, xa, xb};
                            int ind_mesh2[3] = {pair_index[j], xa, xb};
                            atomicAdd(&stress[ind_mesh1], f[xa]*tmp[xb]);
                            atomicAdd(&stress[ind_mesh2], f[xa]*tmp[xb]);
                          }  
                        }
                      }
                      }
                    }
                  }
                }
              }
              ''',
              name='calc_stress')(data, address_list, lmax, sigma, cutoff, L, epsilon, pos, rN, stress)
        stress = cp.asnumpy(stress.reshape(rN, rN, rN, 3, 3))
        local_sigmav += stress/dV
        del stress


    # Time average
    local_mass /= t_ratio
    local_g /= t_ratio
    local_sigma0k /= t_ratio
    local_sigmav /= t_ratio

    for obj in [local_mass, local_g, local_sigma0k, local_sigmav]:
        obj = cp.asnumpy(obj)
   
    local_v = local_g/local_mass[:,:,:,np.newaxis] 
    meshgrid = np.array(rax)
    with h5py.File(outf, 'a') as f:
        f['localMass'][it] = local_mass 
        f['localStress_v'][it] = local_sigma0k - local_g[:,:,:,:,np.newaxis]*local_v[:,:,:,np.newaxis,:]
        f['localStress_f'][it] = local_sigmav
        f['localMomentum'][it] = local_g
        if t == 0:
            f['MeshGrid'] = meshgrid 
        f['time'][it] = time_*timestep *10**-3 / femtosecond

print('Done!')
