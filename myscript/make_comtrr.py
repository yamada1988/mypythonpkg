from pytrr import GroTrrReader
import pytrr
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
trrf0 = GroTrrReader(fname)
with trrf0 as trrfile:
    next(trrfile)
    next(trrfile)
    header_ = trrf0.header
    print(header_)
    endi = header_['endian']
    doub = header_['double']
    tsize = float(header_['time'])
    step = float(header_['step'])
print(endi, doub, tsize, step)

trrf = GroTrrReader(fname)
with trrf as trrfile:
    for data in trrfile:
        data = trrfile.get_data()
        box_size = data['box'][0][0]
        tN += 1

print('Total time:', tN)
print('Read all frame time:', time.time() - t0)

R = make_vec(tN)
V = make_vec(tN) 

t = time.time()
print("dump trrfile to python picklefile (COM information).")
with GroTrrReader(fname) as trrfile:
    it = 0
    for frame in tqdm(trrfile, total=tN):
        data = trrfile.get_data()
        R[it] = [[np.dot(data['x'][3*i:3*i+3],m)*Minv] for i in range(N)]
        V[it] = [[np.dot(data['v'][3*i:3*i+3],m)*Minv] for i in range(N)]
#        print(R[it], V[it])
        it += 1

print('Read trrfile time:', time.time() - t)

t = time.time()
R = np.array(R)
V = np.array(V)
cominfo = {'x':R, 'v':V, 'L': box_size}
ofname = sysname + '_{0:03d}'.format(T) + '.pickle'
with open(ofname, mode='wb') as f:
    pickle.dump(cominfo, f)
print('dump picklefile time:', time.time() - t)

t = time.time()
with open(ofname, mode='rb') as f:
    d = pickle.load(f)
print(d)
print('load picklefile time:', time.time() - t)
