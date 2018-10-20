import h5py
import sys
import math
import numpy as np
import pickle
import time
from tqdm import tqdm

#
# Read data from Gromacs trrfile and write center-of-mass information to picklefile.
# This script requires h5py(https://github.com/andersle/pytrr) and tqdm(https://github.com/tqdm/tqdm).
# Install h5py using "pip insall h5py", "pip install tqdm".
#

def make_vec(tnum):
    vec = ['' for i in range(tN)]
    return vec

args = sys.argv
fname = args[1]
T = int(fname.split('_')[1].split('.')[0])
sysname = fname.split('_')[0]

N = 3000
# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)
Minv = 1/(mO + mH + mH) #(u)
m = np.array([mO, mH, mH])*Minv

tN = 0
t0 = time.time()

f = fname
infh = h5py.File(f, 'r')
tN = len(infh['time'])
box_size = infh['cell_lengths'][0][0]
print('Total time:', tN)
print('Read all frame time:', time.time() - t0)

R = make_vec(tN)
V = make_vec(tN) 

tN = int(tN/2)
t = time.time()
print("dump trrfile to python picklefile (COM information).")

cominfo = {'x':[[] for it in range(tN)], 'v':[[]for it in range(tN)], 'L': box_size}
    
for it in tqdm(range(tN)):
    data_x = infh['coordinates'][it]    
    data_v = infh['velocities'][it]
    R = np.array([np.dot(data_x[4*i:4*i+3].T,m) for i in range(N)])
    V = np.array([np.dot(data_v[4*i:4*i+3].T,m) for i in range(N)])

    cominfo['x'][it] = np.around(R, decimals=4) 
    cominfo['v'][it] = np.around(V, decimals=4)
ofname = sysname + '_{0:03d}'.format(T) + '.pickle'
with open(ofname, mode='wb') as f:
    pickle.dump(cominfo, f)
print('dump picklefile time:', time.time() - t)
