import numpy as np
from scipy.sparse import csr_matrix
from numba import jit
import multiprocessing
import itertools
import time

N = 10000
rN = 12
L =  12.0
dx = L/rN
cutoff = 1.0
n_mesh = rN**3
basic_cell = np.arange(rN)
data = np.arange(n_mesh)
pair_list = np.zeros((rN**3,8),dtype=int)

def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)),
                      shape=(data.max() + 1, data.size))

def get_indices(ind):
    return M[ind].data 


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


@jit('f8[:,:,:](f8[:,:], f8[:,:], i8, i8, f8, f8, f8[:,:,:])', nopython=True)
def calc_dist(vec1, vec2, i1, i2, L, cutoff, array1):
    N = vec1.shape[0]
    M = vec2.shape[0]
    l = 0
    for i in range(N):
        for j in range(M):
            d = 0.0e0
            for k in range(3):
                tmp = vec1[i,k] - vec2[j,k]
                d += np.power((tmp-np.round(tmp/L)*L),2)
            if d < cutoff:
                array1[i1, i2, l] = np.sqrt(d)
            l += 1
 
    return array1


meshfilter = np.array([1,rN,rN**2])
pos = np.random.rand(N,3)*L

t0 = time.time()
p_address = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

# get
# benchmark data: 32x32x32 grid, 1000000 particles
# numpy(single): 31.76 (s)
# cupy(single):   9.35 (s)
# numpy(36process): 1.83 (s)
# use scipy.sparse + list comp : 0.92 (s)
# use scipy.sparse + mp.Pool : 0.34 (s)
t1 = time.time()

M = compute_M(p_address)

p = multiprocessing.Pool()
address_list = p.map(get_indices, data)
p.close()
print('create particle list:',time.time()-t1)

pairs = np.zeros((n_mesh, 8), dtype=np.int64)
t2 = time.time()
[get_pair_list_numba(ind, rN, pairs) for ind in data]
print('create pair_list:',time.time()-t2)

t3 = time.time()
D = np.zeros((n_mesh, 8, 200), dtype=np.float64)
for mesh_ind in data:
    X = pos[address_list[mesh_ind]]
    for ipair, pair in enumerate(pairs[mesh_ind]):
        Y = pos[address_list[pair]]
        #    print("="*28)
        D = calc_dist(X, Y, mesh_ind, ipair, L, cutoff, D)
        #if mesh_ind > 170:
        #    print("-"*32)

print('calc distance:', time.time()-t3)
#for id_, d_ in enumerate(D):
#    print("=== id:",id_)
#    print(address_list[id_])
#    print(pairs[id_])
#    print(d_)
#    print(d_[np.nonzero(d_)])



print('total:', time.time()-t0)
import sys
sys.exit()
#i_ = 0
t3 = time.time()
p = multiprocessing.Pool()
rij = p.map(calc_dist, data)
print('calc distance:', time.time()-t3)
#print(i_/(N**2*0.50))
#print(rij)
