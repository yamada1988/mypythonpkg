import numpy as np
import cupy as cp
import multiprocessing
import time

num_process = 72
N = 1000000
rN = 32
L = 32.0
dx = L/rN

def get_index(ind):
    return np.where(meshids[:] == ind)[0]
#    return cp.where(meshids[:] == ind)[0]   


#meshfilter = cp.array([1,rN,rN**2])
meshfilter = np.array([1,rN,rN**2])

#pos = cp.random.rand(N,3)*L
pos = np.random.rand(N,3)*L

t0 = time.time()
# meshid
meshids = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
#meshids = cp.sum((pos/dx).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

# get
# benchmark data: 32x32x32 grid, 1000000 particles
# numpy(single): 31.76 (s)
# cupy(single):   9.35 (s)
# numpy(36process): 1.83 (s)
t1 = time.time()

# single process:
#data = np.arange(rN**3)
#plist = [get_index(meshids, i) for i in data]

# multi process:
data = np.arange(rN**3)
p = multiprocessing.Pool()
plist = p.map(get_index, data)
print('create particle list:',time.time()-t1)
#print(plist)
