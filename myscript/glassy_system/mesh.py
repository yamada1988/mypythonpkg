import numpy as np
#import cupy as cp
import multiprocessing
import itertools
import time

N = 1000
rN = 4
L =  4.0
dx = L/rN
n_mesh = rN**3
basic_cell = np.arange(rN)
data = np.arange(n_mesh)
pair_list = np.zeros((rN**3,8),dtype=int)
#basic_cell = cp.arange(rN)
#data = cp.arange(n_mesh)
#pair_list = cp.zeros((rN**3,8),dtype=int)

def get_index(ind):
    return np.where(p_address[:] == ind)[0]

def get_pair_list(ind):
    iz = int(ind/(rN**2))
    iy = int((ind-rN**2*iz)/rN)
    ix = (ind-rN**2*iz-rN*iy)
    return np.array([np.roll(basic_cell,ix+dx)[0]+rN*np.roll(basic_cell,iy+dy)[0]+rN**2*np.roll(basic_cell,iz+dz)[0] for dx in range(0,2) for dy in range(0,2) for dz in range(0,2)])

def calc_dist(ind):
    my = address_list[ind]
    return np.array([np.array([np.sqrt(np.sum((pos[pair[0]]-pos[pair[1]])**2)) for ip,pair in enumerate(itertools.product(my,address_list[pl]))]) for pl in pair_list[ind]])

#meshfilter = cp.array([1,rN,rN**2])
meshfilter = np.array([1,rN,rN**2])

#pos = cp.random.rand(N,3)*L
pos = np.random.rand(N,3)*L

t0 = time.time()
# p_address
p_address = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
#p_addres = cp.sum((pos/dx).astype(int)*meshfilter,axis=1)
print('create meshid:',time.time()-t0)

# get
# benchmark data: 32x32x32 grid, 1000000 particles
# numpy(single): 31.76 (s)
# cupy(single):   9.35 (s)
# numpy(36process): 1.83 (s)
t1 = time.time()

# single process:
#plist = [get_index(i) for i in data]

# multi process:
p = multiprocessing.Pool()
address_list = p.map(get_index, data)
print('create particle list:',time.time()-t1)

#for i,a in enumerate((pos/dx).astype(int)):
#    print(i,a)
#for i,p in enumerate(address_list):
#    print(i,p)

t2 = time.time()
p = multiprocessing.Pool()
pair_list = p.map(get_pair_list, data)
print('create pair_list:',time.time()-t2)
#print(pair_list)


#i_ = 0
t3 = time.time()
p = multiprocessing.Pool()
rij = p.map(calc_dist, data)
print('calc distance:', time.time()-t3)
#print(i_/(N**2*0.50))
print(rij)
