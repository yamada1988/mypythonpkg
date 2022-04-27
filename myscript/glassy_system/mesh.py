import numpy as np
from scipy.sparse import csr_matrix
from numba import jit
import multiprocessing
import itertools
import time

N = 160000
rN = 32
L =  32.0
dx = L/rN
cutoff = 0.850
n_mesh = rN**3
N1 = int(N/2.0)
s1 = 1.0e0
s2 = 1.20e0
sigma1 = s1 * np.ones(N1, dtype=np.float64)
sigma2 = s2 * np.ones(N1, dtype=np.float64)
sigma = np.hstack((sigma1, sigma2))
basic_cell = np.arange(rN)
data = np.arange(n_mesh)
pair_list = np.zeros((rN**3,8),dtype=int)

def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)),
                      shape=(data.max() + 1, data.size))

def get_indices(ind):
    return M[ind].data.tolist()


@jit('i8[:,:](i8, i8, i8[:, :])', nopython=True)
def get_pair_list_numba(ind, rN, pairs):
    iz = int(ind/(rN**2))
    iy = int((ind-rN**2*iz)/rN)
    ix = int(ind-(rN**2*iz+rN*iy))
    j = 0
    for dz in range(2):
        for dy in range(2):
            for dx in range(2):
                ix_ = (ix + dx) % rN
                iy_ = (iy + dy) % rN
                iz_ = (iz + dz) % rN
                pair = ix_ + iy_*rN + iz_*rN**2
                pairs[ind, j] = pair
                j += 1
    return pairs


@jit('f8[:,:,:](f8[:,:], i8[:,:], i8[:,:], i8, i8, f8, f8, f8[:], f8[:,:,:])', nopython=True)
def calc_stress(pos, address_list, pair_index, n_mesh, N1, L, cutoff, sigma, array1):
    ep = 1.0e0
    c2 = cutoff**2
    for mesh_ind in range(n_mesh):
        i1s = address_list[mesh_ind]
        i_ = 1
        for i1 in i1s:
            if i1 < 0:
                continue
            X = pos[i1]
            for ipair in range(8):
                if ipair == 0:
                    i2s = address_list[pair_index[mesh_ind][ipair]][i_:] # avoid double count in same grid
                else:
                    i2s = address_list[pair_index[mesh_ind][ipair]]

                for i2 in i2s:
                    if i2 < 0:
                        continue
                    Y = pos[i2]
                    d = 0.0e0
                    tmp = np.zeros(3, dtype=np.float64)
                    for k in range(3):
                        tmp[k] = X[k] - Y[k]
                        tmp[k] -= np.round(tmp[k]/L)*L
                        d += np.power(tmp[k],2)
           
                    if d < c2:
                        r = np.sqrt(d)
                        rinv = np.power(0.50e0*(sigma[i1]+sigma[i2])/r, 6)
                        f = 4.0e0*ep*6.0e0*rinv*(-2.0e0*rinv+1.0e0)*0.50e0*tmp/r 
                        temp_s = 0.50e0*np.outer(f, tmp)
                        array1[mesh_ind] += temp_s
                        array1[pair_index[mesh_ind][ipair]] += temp_s
            i_ += 1
    return array1


meshfilter = np.array([1,rN,rN**2])
pos = np.random.rand(N,3)*L

t0 = time.time()
p_address = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

t1 = time.time()

M = compute_M(p_address)
p = multiprocessing.Pool()
X = p.map(get_indices, data)
p.close()
## Max size
l = max(map(len, X))
## Padding
address_list = np.array(map(lambda x: x + [-1]*(l-len(x)), X), dtype=np.int64)
print('create particle list:',time.time()-t1)

pair_index = np.zeros((n_mesh, 8), dtype=np.int64)
t2 = time.time()
[get_pair_list_numba(ind, rN, pair_index) for ind in data]
print('create pair_list:',time.time()-t2)

t3 = time.time()
stress_tensor = calc_stress(pos, address_list, pair_index, n_mesh, N1, L, cutoff, sigma, np.zeros((n_mesh, 3, 3), dtype=np.float64))
print(stress_tensor)
print('calc stress:', time.time()-t3)

print('total:', time.time()-t0)

import sys
sys.exit()
