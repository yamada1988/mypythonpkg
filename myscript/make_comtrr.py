from pytrr import GroTrrReader
import sys
import math
import numpy as np
import pickle
import time
from tqdm import tqdm
#
# Read data from Gromacs trrfile and write center-of-mass information to picklefile.
# This script requires pytrr(https://github.com/andersle/pytrr) and tqdm(https://github.com/tqdm/tqdm).
# Install pytrr using "pip insall pytrr", "pip install tqdm".
#

def make_vec(tnum, inum):
    vec = [['' for i in range(tN)] for i in range(N)]
    return vec

args = sys.argv
fname = args[1]
T = int(fname.split('_')[1].split('.')[0])
sysname = fname.split('_')[0]

N = 3000
# Mass Parameters
mO = 16.00 * 1.661e-27 #(kg)
mH =  1.08 * 1.661e-27 #(kg)
M = mO + mH + mH #(u)

tN = 0
t0 = time.time()
with GroTrrReader(fname) as trrfile:
    for frame in trrfile:
        frame_data = trrfile.get_data()
        tN += 1
print('Total time;', tN)
print('Read all frame time:', time.time() - t0)

R = make_vec(tN, N)
V = make_vec(tN, N) 

t = time.time()
print("dump trrfile to python picklefile (COM information).")
with GroTrrReader(fname) as trrfile:
    it = 0
    for frame in tqdm(trrfile, total=tN):
        if it >= tN:
            break
        data = trrfile.get_data()
        for i in range(N):
            iO  = 3*i 
            iH1 = 3*i + 1
            iH2 = 3*i + 2

# Calculate single H2O molecule's center of mass (R)
            R[i][it] =(mO*data['x'][iO]+mH*data['x'][iH1]+mH*data['x'][iH2])/M

# Calculate single H2O molecule's velocity (V)
            V[i][it] =(mO*data['v'][iO]+mH*data['v'][iH1]+mH*data['v'][iH2])/M
        it += 1

print('Read trrfile time:', time.time() - t)

cominfo = {'x':R, 'v': V}
ofname = sysname + '_{0:03d}'.format(T) + '.pickle'
with open(ofname, mode='wb') as f:
    pickle.dump(cominfo, f)
