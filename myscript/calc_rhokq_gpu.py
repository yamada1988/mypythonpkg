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
# python script/calc_rhokq.py MD/sample0001_01.pickle 1 5 0 5
#

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nens = fname.split('_')[0][-4:]
dirname = fname.split('/')[0]
sysname = fname.split('/')[1].split('_')[0]
nk_min = int(args[2])
nk_max = int(args[3])
nq_min = int(args[4])
nq_max = int(args[5])

pi = math.pi

t = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)

tN = len(d['x'])
talpha = 0.010e0
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
kN = nk_max - nk_min + 1

t0 = time.time()

K = np.array([1/math.sqrt(3.0e0) * np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(nk_min, nk_max+1)])
Q = np.array([1/math.sqrt(3.0e0) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(nq_min, nq_max+1)])
K = np.around(K, decimals=4)
Q = np.around(Q, decimals=4)
K = cp.asarray(K)
Q = cp.asarray(Q)

t1 = time.time()
bnum = 200
tN_b = int(tN/bnum)
rho = np.array([[(0.0e0+0.0e0j) for q in Q] for k in K])
theta = np.array([[np.zeros((tN_b, N)) for q in Q] for k in K])
rho = cp.asarray(rho)
theta = cp.asarray(theta)
print('rho:')
print(rho.shape)

# Calculate rho(k,q,t)
R = d['x'][:tN]
R = cp.asarray(R)

print(R.shape)
print('theta(r,k)')
for it in range(tN_b):
    r = R[it*bnum:(it+1)*bnum]
    for ik,k in enumerate(K):
        for iq,q in enumerate(Q):
            a = cp.dot(r, k)
            theta[ik][iq][it] = a[it]
del R
del a

V = d['v'][:][:tN]
V = cp.asarray(V)

print('theta(v,q)')
for ib,it in enumerate(range(tN_b)):
    v = V[ib*bnum:(ib+1)*bnum]
    for ik,k in enumerate(K):
        for iq,q in enumerate(Q):
            a = cp.dot(v, q)
            theta[ik][iq][it] += a[it]
del V


for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        print('(ik,iq)=({0:d},{1:d})'.format(ik,iq))
        a = cp.sum(cp.exp(-1.0e0j*theta),axis=3)
        b = cp.sum(a, axis=2) /  tN_b
        print('b:', b.shape)
        print(b[ik][iq])
        rho[ik][iq] = b[ik][iq]
        print('<rho(k,q,t)>:')
        print(rho[ik][iq])
del a
del b

print("Calculate rho(k,q,t) time:", time.time()-t1)
rho = cp.asnumpy(rho)
t3 = time.time()
for ik,nk in enumerate(range(nk_min, nk_max+1)):
    ofname = 'DAT/rho'+nens+'nk{0:02d}_{1:d}.dat'.format(nk,int(T))
    with open(ofname, 'wt') as f:
        f.write('# rho[k][0] = {0}, {1}\n'.format(rho[ik][0].real, rho[ik][0].imag))
        f.write('# k=[{0[0]},{0[1]},{0[2]}]\n# rho(k,q,t)\n'.format(K[ik]))
    for iq,q in enumerate(Q):
        a = rho[ik][iq].real/rho[ik][0].real
        b = rho[ik][iq].imag/rho[ik][0].imag
        with open(ofname, 'a+') as f:
            f.write('{0}\t{1}\t{2}\n'.format(q[0], a, b))

print("write file time:", time.time()-t3)
print("total time:", time.time()-t0)
