import numpy as np
import time

N = 1000000
rN = 16
L = 32.0
dx = L/rN

def get_index(array, ind):
    return np.where(array[:] == i)[0]
    

meshfilter = np.array([1,rN,rN**2])

pos = np.random.rand(N,3)*L

t0 = time.time()
# meshid
meshids = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

# get
t1 = time.time()
plist = [get_index(meshids, i) for i in range(rN**3)]
print('create particle list:',time.time()-t1)
print(plist)
