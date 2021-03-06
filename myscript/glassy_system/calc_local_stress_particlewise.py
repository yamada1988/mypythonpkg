import cupy as cp
import numpy as np
import time 

N = 160000
N1 = N/2
L =  32.0
cutoff = 0.850
rN = int(L/cutoff)
n_mesh = rN**3
data = cp.arange(N, dtype=np.int32)
s1 = 0.340e0
s2 = 0.408e0
epsilon = 0.994
sigma1 = s1*cp.ones(N1)
sigma2 = s2*cp.ones(N1)
sigma = cp.concatenate((sigma1, sigma2))
sigma = cp.asarray(sigma, dtype=np.float32)

pos = L*cp.random.rand(N,3, dtype=np.float32)
pos = cp.asarray(pos, dtype=np.float32)
r0 = cp.asarray(pos[:,0]/cutoff, dtype=np.int32)
r1 = cp.asarray(pos[:,1]/cutoff, dtype=np.int32)
r2 = cp.asarray(pos[:,2]/cutoff, dtype=np.int32)
stress = cp.zeros((rN, rN, rN, 3, 3), dtype=np.float32)

t0 = time.time()
lmax = 2000
int_list = -1*cp.ones((N, lmax), dtype=np.int32)
update_list = cp.ElementwiseKernel(
    in_params='int32 data, raw float32 pos, float32 cutoff, float32 L', 
    out_params='raw int32 int_list',
    operation=\
    '''
    int l = 0;
    for (int j = i+1; j < _ind.size(); j++){
         double tmp[3] = {0.0};
         double diff = 0.0;

         for (int k = 0; k<3; k++){ 
           tmp[k] = pos[3*j+k] - pos[3*i+k];
           if (tmp[k] < -0.50*L) {
           tmp[k] += L;
           } else if (tmp[k] > 0.50*L) {
           tmp[k] -= L;
           } 
           diff += powf(tmp[k],2);
         }
         double r = sqrt(diff);

         if ( r > 5.0 *cutoff) {
           continue;
         } else{
                 int ind[] = {i, l};
                 int_list[ind] = j;
                 l += 1;
           }
    }
    ''',
    name='update_list')(data, pos, cutoff, L, int_list)

print('update list:', time.time()-t0)

t1 = time.time()
calc_stress = cp.ElementwiseKernel(
    in_params='int32 data, raw int32 r0, raw int32 r1, raw int32 r2, raw float32 sigma, float32 cutoff, float32 L, float32 epsilon, raw float32 pos, raw int32 int_list, int32 lmax', 
    out_params='raw float32 stress',
    operation=\
    '''
    for (int j = 0; j < lmax; j++){
         int jp = int_list[lmax*i+j];
         if ( jp<0){
         break;
         }
         else
         {
         double tmp[3] = {0.0};
         double diff = 0.0;

         for (int k = 0; k<3; k++){ 
           tmp[k] = pos[3*jp+k] - pos[3*i+k];
           tmp[k] -= L*int(tmp[k]/(L/2.0));
           diff += powf(tmp[k],2);
         }
         double r = sqrt(diff);

         if ( r > cutoff ) {
           continue;
         } else{
             double r6 = powf(0.50*(sigma[i]+sigma[j])/r, 6);
             double f0 = -24.0*epsilon*r6*(2.0*r6-1.0)/r;
             for (int l=0; l<3; l++){
               for (int m=0; m<3; m++){
                 int ind_mesh1[5] = {r0[i], r1[i], r2[i], l, m};
                 int ind_mesh2[5] = {r0[j], r1[j], r2[j], l, m};
                 atomicAdd(&stress[ind_mesh1], f0*tmp[l]*tmp[m]);
                 atomicAdd(&stress[ind_mesh2], f0*tmp[l]*tmp[m]);
               }
             }
           }
         }
    }
    ''',
    name='calc_stress')(data, r0, r1, r2, sigma, cutoff, epsilon, L, pos, int_list, lmax,  stress)

print('calc stress:', time.time()-t1)
print('stress:')
stress = cp.asnumpy(stress)
print(stress)
import sys
sys.exit()
