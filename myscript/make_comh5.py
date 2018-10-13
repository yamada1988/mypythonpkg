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

def calc_RV(list_data, i, it, m):
    # Calculate single H2O molecule's center of mass (R)
    #R[i][it] = (mO*list_data['x'][iO]+mH*list_data['x'][iH1]+mH*list_data['x'][iH2])/M
    r = np.dot(list_data['x'][3*i:3*i+3],m)/M

    # Calculate single H2O molecule's velocity (V)
    #V[i][it] = (mO*list_data['v'][iO]+mH*list_data['v'][iH1]+mH*list_data['v'][iH2])/M
    v = np.dot(list_data['v'][3*i:3*i+3],m)/M
    return r, v    

args = sys.argv
fname = args[1]
T = int(fname.split('_')[1].split('.')[0])
sysname = fname.split('_')[0]

N = 3000
# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)
Minv = 1/(mO + mH + mH) #(u)
m = np.array([mO, mH, mH])

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

t = time.time()
print("dump trrfile to python picklefile (COM information).")
for it in tqdm(range(tN)):
    data_x = infh['coordinates'][it]    
    data_v = infh['velocities'][it]
    R[it] = [np.dot(data_x[3*i:3*i+3],m)*Minv for i in range(N)]
    V[it] = [np.dot(data_v[3*i:3*i+3],m)*Minv for i in range(N)]
print('Read trrfile time:', time.time() - t)

t = time.time()
R = np.array(R)
V = np.array(V)
cominfo = {'x':R, 'v':V, 'L': box_size}
ofname = sysname + '_{0:03d}'.format(T) + '.pickle'
with open(ofname, mode='wb') as f:
    pickle.dump(cominfo, f, protocol=4)
print('dump picklefile time:', time.time() - t)
