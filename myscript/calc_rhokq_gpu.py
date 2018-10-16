import sys
import math
import numpy as np
import pickle
import statsmodels.api as sm
import cupy as cp
import tqdm
import time


#
# Calculate rho(k,q,t), <rho(k,q)>, and <drho(k,q,t)drho(-k,-q,0)> from compickle.
# This script requires cupy(https://github.com/cupy/cupy) 
# and statsmodels(https://github.com/statsmodels/statsmodels).
# Install cupy using "conda install cupy".
# Install statsmodels using "using install statsmodels".
#

# Sample:
# python script/calc_rhokq.py MD/sample0001_01.pickle 1 0 5
#

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nens = fname.split('_')[0][-4:]
dirname = fname.split('/')[0]
sysname = fname.split('/')[1].split('_')[0]
nk = int(args[2])
nq_min = int(args[3])
nq_max = int(args[4])

pi = math.pi

t = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)

tN = len(d['x'])
talpha = 1.0e0
tmax = int(tN*talpha)
tN = tmax
box_size = d['L']
R = d['x']

print(tN, box_size)
print('load picklefile time:', time.time() - t)


# Phisical Parameters
N =  len(R[0])
dt = 1.0e0 #(ps)
L = box_size #(nm) 
mO = 16.00e-27 #(kg)
mH =  1.008e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)

print('R:')
print(len(R))
del R
# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 1.0e0*2.0e0*pi / L # (nm^-1)
vth = math.sqrt(3.0e0*kB*T*Minv) /1000.0e0 # (nm/fs)^-1 = (km/s)^-1
alpha = 0.0250
dq = alpha * 2.0e0*pi / vth
q0 = 0.0e0
qN = nq_max - nq_min + 1

t0 = time.time()

K = np.array([1/math.sqrt(3.0e0) * np.array([k0+nk*dk, k0+nk*dk, k0+nk*dk]) for i in range(nq_min, nq_max+1)])
Q = np.array([1/math.sqrt(3.0e0) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(nq_min, nq_max+1)])
K = np.around(K, decimals=4)
Q = np.around(Q, decimals=4)
K = cp.asarray(K)
Q = cp.asarray(Q)

t1 = time.time()
bnum = 100
tN_b = int(tN/bnum)
print('tN_b:',tN_b)
rho = np.array([(0.0e0+0.0e0j) for q in Q])
rho = cp.asarray(rho)
print('rho:')
print(rho.shape)

# Calculate rho(k,q,t)
R = d['x'][:tN]
V = d['v'][:tN]
print(K.shape, Q.shape)
print('theta(r,k)')

c = np.array([0.0e0+0.0e0j for i in range(qN)])
b = cp.asarray(c)
for it in range(tN_b):
    r = cp.asarray(R[it*bnum:(it+1)*bnum])
    v = cp.asarray(V[it*bnum:(it+1)*bnum])
    a = cp.dot(r, K.T) + cp.dot(v, Q.T)
    b += cp.sum(cp.exp(-1.0e0j*a),axis=(0,1))
b /= float(tN_b)
del R
del V
for iq,q in enumerate(Q):
    print('b:', b.shape)
    print(b[iq])
    rho[iq] = b[iq]
    print('<rho(k,q,t)>:')
    print(rho[iq])
del b
del c

print("Calculate rho(k,q,t) time:", time.time()-t1)
rho = cp.asnumpy(rho)
t3 = time.time()
ofname = 'DAT/rho'+nens+'nk{0:02d}_{1:d}.dat'.format(nk,int(T))
with open(ofname, 'wt') as f:
    f.write('# rho[k][0] = {0}, {1}\n'.format(rho[0].real, rho[0].imag))
    f.write('# k=[{0[0]},{0[1]},{0[2]}]\n# rho(k,q,t)\n'.format(K[0]))
for iq,q in enumerate(Q):
    a = rho[iq].real/rho[0].real
    b = rho[iq].imag/rho[0].imag
    with open(ofname, 'a+') as f:
        f.write('{0}\t{1}\t{2}\n'.format(q[0], a, b))

print("write file time:", time.time()-t3)
print("total time:", time.time()-t0)
