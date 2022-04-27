from pytrr import GroTrrReader
import sys
import math
import numpy as np
import pickle
import statsmodels.api as sm
import tqdm
import time


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
R = d['x'][:3000]
V = d['v'][:3000]
box_size = d['L']
print(tN, box_size)
print('load picklefile time:', time.time() - t)

# Phisical Parameters
N =  len(R[0])
dt = 1.0e0 #(ps)
L = box_size #(nm) 
mO = 16.00*1.66e-27 #(kg)
mH =  1.008*1.66e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 1.0e0*2.0e0*pi / L # (nm^-1)
vth = math.sqrt(kB*T*Minv) /1000.0e0 # (nm/fs)^-1 = (km/s)^-1
alpha = 0.0250
dq = alpha * 2.0e0*pi / vth
q0 = 0.0e0
qN = nq_max - nq_min + 1
kN = nk_max - nk_min + 1

t0 = time.time()

#K = np.array([1/math.sqrt(3) *np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(nk_min, nk_max+1)])
#Q = np.array([1/math.sqrt(3) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(nq_min, nq_max+1)])
K = np.array([[k0 + dk*ix, k0 + dk*iy, k0 + dk*iz] for ix in range(kN) for iy in range(kN) for iz in range(kN)])
Q = np.array([[q0 + dq*ix, q0 + dq*iy, q0 + dq*iz] for ix in range(qN) for iy in range(qN) for iz in range(qN)])

# Zip R[it][i] and V[it][i]
RV = np.dstack((R,V))
for k in K:
    for i,ki in enumerate(k):
        ki = round(ki, 6)
        k[i] = ki
for q in Q:
    for i,qi in enumerate(q):
        qi = round(qi, 6)
        q[i] = qi
# Zip K[ik] and Q[iq]
KQ = np.array([ (k,q) for q in Q for k in K])
print(KQ.shape)
KQ = KQ.reshape(kN**3*qN**3,6)

t1 = time.time()
# Calculate rho(k,q,t)
theta = np.dot(RV,KQ.T)
rho = np.sum(np.exp(theta*-1j), axis=1)
print("Calculate rho(k,q,t) time:", time.time()-t1)

print(rho.shape)
# Calculate <rho(k,q)>
rho = np.sum(rho, axis=0) / float(tN)

print('rho:', rho.shape)
t3 = time.time()
ofname = 'DAT/rho'+nens+'kq_{0:d}.dat'.format(int(T))
with open(ofname, 'wt') as f:
    f.write('# kx\tky\tkz\tqx\tqy\tqz\trho(k,q,t)\n')
with open(ofname, 'a+') as f:
    for ikq,kq in enumerate(KQ):
        f.write('{0:5.3f}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\t{4:5.3f}\t{5:5.3f}\t{6:6.4f}\t{7:6.4f}\n'.format(kq[0],kq[1],kq[2],kq[3],kq[4],kq[5], rho[ikq].real, rho[ikq].imag))

print("write file time:", time.time()-t3)
print("total time:", time.time()-t0)
