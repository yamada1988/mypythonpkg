from pytrr import GroTrrReader
import sys
import math
import numpy as np
import pickle
import statsmodels.api as sm
import time


#
# Time measuring version
#
# time.sh: 
##!/bin/bash

#
# Calculate rho(k,q,t), <rho(k,q)>, and <drho(k,q,t)drho(-k,-q,0)> from Gromacs trrfile.
# This script requires pytrr(https://github.com/andersle/pytrr) 
# and statsmodels(https://github.com/statsmodels/statsmodels).
# Install pytrr using "pip insall pytrr".
# Install statsmodels using "using install statsmodels".
#

# Sample:
# python script/calc_rhokq.py MD/sample0001_01.pickle 1 5 0 5
#

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nens = fname.split('_')[0][-4:]
nk_min = int(args[2])
nk_max = int(args[3])
nq_min = int(args[4])
nq_max = int(args[5])

pi = math.pi

t = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)

tN = len(d['x'])
R = d['x']
V = d['v']
box_size = d['L']
print(tN, box_size)
print('load picklefile time:', time.time() - t)

# Phisical Parameters
N =  len(R[0])
dt = 0.010e0 #(ps)
L = box_size #(nm) 
mO = 16.00e-27 #(kg)
mH =  1.008e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 1.0*2.0e0*pi / L # (nm^-1)
vth = math.sqrt(3.0e0*kB*T*Minv) # (nm/fs)^-1 = (km/s)^-1
alpha = 50.0
dq = alpha * 2.0e0*pi / vth
q0 = 0.0e0
qN = nq_max - nq_min + 1
kN = nk_max - nk_min + 1

t0 = time.time()

K = np.array([1/math.sqrt(3) *np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(nk_min, nk_max+1)])
Q = np.array([1/math.sqrt(3) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(nq_min, nq_max+1)])
for k in K:
    for i,ki in enumerate(k):
        ki = round(ki, 6)
        k[i] = ki
for q in Q:
    for i,qi in enumerate(q):
        qi = round(qi, 6)
        q[i] = qi

t1 = time.time()
rho = np.array([[[(0+1j) for it in range(tN)] for q in Q] for k in K])
print('rho:')
print(rho.shape)
# Calculate rho(k,q,t)
for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        theta = np.dot(R,k) + np.dot(V,q)
        rho[ik][iq] = np.sum(np.exp(theta*-1j), axis=1)
print("Calculate rho(k,q,t) time:", time.time()-t1)

t2 = time.time()
# Calculate C(k,q,t) = <drho(k,q,t)drho(-k,q,0)>
thalf = int(tN*0.8)
C_kqt = [[[0.0 for it in range(tN)] for q in Q] for k in K]
for ik in range(kN):
    for iq in range(qN):
        C_kqt[ik][iq] = sm.tsa.stattools.acf(rho[ik][iq], nlags=thalf)

print("Calculate C(k,q,t) time:", time.time()-t2)

t3 = time.time()
for ik,nk in enumerate(range(nk_min, nk_max+1)):
    for iq,nq in enumerate(range(nq_min, nq_max+1)):
        ofname = 'DAT/'+nens+'nk{0:02d}nq{1:02d}_{2:d}.dat'.format(nk,nq,int(T))
        with open(ofname, 'wt') as f:
            f.write('# k=[{0[0]},{0[1]},{0[2]}]\n# q=[{1[0]},{1[1]},{1[2]}]\n# t(ps) C(k,q,t)\n'.format(K[ik],Q[iq]))

        for t_interval in range(thalf-1):
            t = '{0:5.2f}'.format(t_interval * dt)
            with open(ofname, 'a+') as f:
                f.write('{0:5.2f}\t{1:9.7f}\t{2:9.7f}\n'.format(float(t), C_kqt[ik][iq][t_interval].real, C_kqt[ik][iq][t_interval].imag))

print("write file time:", time.time()-t3)
print("total time:", time.time()-t0)
