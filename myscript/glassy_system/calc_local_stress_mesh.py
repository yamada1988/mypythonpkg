import numpy as np
import cupy as cp
import time

N = 100
L =  3.0
cutoff = 1.0
rN =int(L/cutoff)
n_mesh = rN**3
N1 = int(N/2.0)
s1 = 0.340e0
s2 = 0.408e0
epsilon = 0.994
sigma1 = s1 * np.ones(N1, dtype=np.float64)
sigma2 = s2 * np.ones(N1, dtype=np.float64)
sigma = np.hstack((sigma1, sigma2))
basic_cell = np.arange(rN)
data = np.arange(n_mesh)
pair_list = np.zeros((rN**3,8),dtype=int)

meshfilter = np.array([1,rN,rN**2])
pos = np.random.rand(N,3)*L

t0 = time.time()
p_address = np.sum((pos/cutoff).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

t1 = time.time()

lmax = 12
address_list = -1*cp.ones((n_mesh, lmax), dtype=np.int32)
r0 = cp.asarray(pos[:,0]/cutoff, dtype=cp.int32)
r1 = cp.asarray(pos[:,1]/cutoff, dtype=cp.int32)
r2 = cp.asarray(pos[:,2]/cutoff, dtype=cp.int32)

data = cp.asarray(data, dtype=np.int32)
cp.ElementwiseKernel(
    in_params='int32 data, raw int32 r0, raw int32 r1, raw int32 r2, float32 cutoff, int32 rN, int32 N, int32 lmax', 
    out_params='raw int32 address_list',
    operation=\
    '''
    int l = 0;
    for (int j = 0; j < N; j++){
         if ( i == r0[j] + rN*r1[j] + rN*rN*r2[j] ) {
           int ind[2] = {i, l};
           address_list[ind] = j;
           l += 1;
         } else {
         continue;
         }
    }
    ''',
    name='update_list')(data, r0, r1, r2, cutoff, rN, N, lmax, address_list)

for ia, a in enumerate(address_list):
    print(ia, a)
#import sys
#sys.exit()

print('create particle list:',time.time()-t1)

t2 = time.time()

sigma = cp.asarray(sigma, dtype=np.float32)
pos = cp.asarray(pos, dtype=np.float32)
stress = cp.zeros((n_mesh, 3, 3), dtype=np.float32)
calc_stress = cp.ElementwiseKernel(
    in_params='int32 data, raw int32 address_list, int32 lmax, raw float32 sigma, float32 cutoff, float32 L, float32 epsilon, raw float32 pos, int32 rN', 
    out_params='raw float32 stress',
    operation=\
    '''
    int n_mesh = _ind.size();
    int pair_index[8] = {i, (i+1)%n_mesh, (i+rN)%n_mesh, (i+rN+1)%n_mesh, (i+rN*rN)%n_mesh, (i+1+rN*rN)%n_mesh, (i+rN+rN*rN)%n_mesh, (i+1+rN+rN*rN)%n_mesh};
    for (int p1=0; p1<lmax; p1++){
      int ip = address_list[lmax*i+p1];
      if ( ip < 0 )
      {
      break;
      } 
      else 
      {
      for (int p2=0; p2<lmax; p2++){
        for (int j=0; j<8; j++){
          int jp = address_list[lmax*pair_index[j]+p2];
          if ( jp<0 || ip == jp )
          {
          break;
          }
          else
          {
          double tmp[3] = {0.0};
          double diff = 0.0;
          for (int m = 0; m<3; m++){
            tmp[m] = pos[3*ip+m] - pos[3*jp+m];
            if (tmp[m] < -0.50*L) {
            tmp[m] += L;
            } else if (tmp[m] > 0.50*L) {
            tmp[m] -= L;
            } 
            diff += powf(tmp[m],2); 
          }
          double r = sqrt(diff);
          if ( r > cutoff)
          {
          continue;
          }
          else
          {
            double r6 = powf(0.50*(sigma[ip]+sigma[jp])/r, -6);
            double f0 = -24.0*epsilon*r6*(2*r6-1.0)/r;
            for (int xa=0; xa<3; xa++){
              for (int xb=0; xb<3; xb++){
                int ind_mesh1[3] = {i, xa, xb};
                int ind_mesh2[3] = {pair_index[j], xa, xb};
                atomicAdd(&stress[ind_mesh1], f0*tmp[xa]*tmp[xb]);
                atomicAdd(&stress[ind_mesh2], f0*tmp[xa]*tmp[xb]);
                }  
              }
          }
          }
          }
      }
      }
    }
    ''',
    name='calc_stress')(data, address_list, lmax, sigma, cutoff, epsilon, L, pos, rN, stress)

print('calc stress:', time.time()-t2)
print('total:', time.time()-t0)
#stress = cp.reshape(stress, (rN,rN,rN,3,3))
for is_, s in enumerate(stress):
    print(is_) 
    print(s)

import sys
sys.exit()
