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
kqN = qN**3*kN**3

t0 = time.time()

#K = np.array([1/math.sqrt(3) *np.array([k0+i*dk, k0+i*dk, k0+i*dk]) for i in range(nk_min, nk_max+1)])
#Q = np.array([1/math.sqrt(3) * np.array([q0+i*dq, q0+i*dq, q0+i*dq]) for i in range(nq_min, nq_max+1)])
K = np.array([[k0 + dk*ix, k0 + dk*iy, k0 + dk*iz] for ix in range(kN) for iy in range(kN) for iz in range(kN)])
Q = np.array([[q0 + dq*(ix-qN/2), q0 + dq*(iy-qN/2), q0 + dq*(iz-qN/2)] for ix in range(qN) for iy in range(qN) for iz in range(qN)])
K = np.around(K, decimals=4)
Q = np.around(Q, decimals=4)


t1 = time.time()
bnum = 1
tN_b = int(tN/bnum)
print('tN_b:',tN_b)
# Zip R[it][i] and V[it][i]
RV = np.dstack((R,V))
# Zip K[ik] and Q[iq]
KQ = np.array([ (k,q) for q in Q for k in K])
print(KQ.shape)
KQ = KQ.reshape(kqN,6)

# Set cupyarray
KQ = cp.asarray(KQ)
RV = cp.asarray(RV)

print(RV.shape)

b = np.array([0.0e0+0.0e0j for ikq in range(kN**3*qN**3)])
print('b:', b.shape)
rho = np.array([0+0j for ikq in range(kN**3*qN**3)])
b = cp.asarray(b)
rho = cp.asarray(rho)
for it in range(tN_b):
    print('it:',it)
    rv = cp.asarray(RV[it*bnum:(it+1)*bnum])
    #print(rv.shape,KQ.shape)
    a = cp.dot(rv, KQ.T)
    test = cp.exp(-1.0e0j*a)
    test2 = cp.sum(test,axis=(0,1))
    #print('test2:', test2.shape)
    rho += test2
rho /= float(tN_b)
del RV
del b

print(rho)
print('rho:', rho.shape)
for ikq,kq in enumerate(KQ):
    print(rho[ikq])

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
