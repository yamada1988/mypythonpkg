import sys
import math
import numpy as np
import pickle
from scipy.spatial import distance
from scipy import integrate
import time


#
# This script calculate S(k,q) from compickle.
# Sample:
# python script/calc_Skq.py MD/sample0001_01.pickle 1 5 0 5
#

def func(r,v):
    r_box = int(math.ceil(r/dr))
    v_box = int(math.ceil(v/dv))
    if k == 0.0e0 and not q == 0.0e0:
        dummy = P[r_box][v_box]*(2.0e0*r**2)*(v*np.sin(q*v)/(q))*(2.0e0*math.pi)**2
    elif not k == 0.0e0 and q == 0.0e0:
        dummy = P[r_box][v_box]*(r*np.sin(k*r)/(k))*(2.0e0*v**2)*(2.0e0*math.pi)**2
    elif k == 0.0e0 and q == 0.0e0:
        dummy = P[r_box][v_box]*(2.0e0*r**2)*(2.0e0*v**2)*(2.0e0*math.pi)**2
    else:
        dummy = P[r_box][v_box]*(r*np.sin(k*r)/(k))*(v*np.sin(q*v)/(q))*(2.0e0*math.pi)**2
    return dummy 

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

print('load picklefile:')
t = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)
print('time:', time.time() - t)

tN = len(d['x'])
talpha = 0.40e0
tmax = int(tN*talpha)
tN = tmax
box_size = d['L']
R = d['x']
print(tN, box_size)

# Phisical Parameters
N =  len(R[0])
dt = 1.0e0 #(ps)
L = box_size #(nm) 
mO = 16.00e-27 #(kg)
mH =  1.008e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)

# Wave number Parameters
k0 = 2.0e0*pi / L # (nm^-1)
dk = 1.0e0*2.0e0*pi / L # (nm^-1)
kN = nk_max - nk_min + 1
K = [k0+ik*dk for ik in range(kN)]
vth = math.sqrt(3.0e0*kB*T*Minv) / 1000.0e0 # (nm/ps)^-1 = (m/s)^-1
alpha = 0.0050
dq = alpha * 2.0e0*pi / vth
q0 = 0.0e0
qN = nq_max - nq_min + 1
Q = [q0 + iq*dq for iq in range(qN)]

# Space-Velocity Parameters
r_min = 0.0e0
r_max = float(L/2.0e0)
rN = 100
dr = (r_max-r_min)/ float(rN)
r_ = np.array([r_min + ir*dr for ir in range(rN+1)])
v_min = 0.0e0
v_max = 12.0e0
vN = 50
dv    = (v_max-v_min)/ float(vN)
v_ = np.array([v_min + iv*dv for iv in range(vN+1)])

bnum = 1
tN_b = int(tN/bnum)
totalN = int(N*(N-1)/2)
print('tN_b:',tN_b)

print(r_max)
# Calculate rho(k,q,t)
R = d['x']
V = d['v']
P = [[0.0e0 for iv in range(vN+1)] for j in range(rN+1)]
for it in range(tN_b):
    print('it:', it)
    rvec = np.array(R[it])
    vvec = np.array(V[it])
    # Cannot use distance.pdist due to PBC
    xr = rvec[np.newaxis, :] - rvec[:, np.newaxis]
    xr -= L * np.trunc(xr/L) # unset PBC
    rij = np.sqrt(np.sum(xr**2, axis=2))
    r = rij[np.triu_indices(N, k = 1)]
    vij = distance.pdist(vvec, 'euclidean')
    v = vij.flatten()

    for ir,r_dummy in enumerate(r):
        v_dummy = v[ir]
        r_box = int(math.ceil(r_dummy/dr))
        v_box = int(math.ceil(v_dummy/dv))
        if r_box <= rN:
            P[r_box][v_box] += 1
print('sanity check:')
P = np.array(P)
count_total = np.sum(P)
print('count_total:{0:d} totalN:{1:d}'.format(int(count_total), int(totalN*tN)))
# normalization
P /= float(count_total)
dVr = np.zeros(rN)
dVv = np.zeros(vN)
dVr = 4.0e0*pi*(r_+dr/2.0e0)**2*dr
dVv = 4.0e0*pi*(v_+dv/2.0e0)**2*dv
for i in range(rN):
    for j in range(vN):
        P[i][j] /= (dVr[i]*dVv[j])

# Integration
S = np.zeros((kN, qN))
for ik,k in enumerate(K):
    for iq,q in enumerate(Q):
        S[ik][iq], _ = integrate.nquad(func, [[0.0, r_max], [0.0, v_max]])  
        print(S[ik][iq])

    ofname = 'DAT/Skq'+nens+'nk{0:02d}_{1:d}.dat'.format(ik,int(T))
    with open(ofname, 'wt') as f:
        f.write('# S[k][0] = {0}\n'.format(S[ik][0]))
        f.write('# k={0:8.5f}\n# S(k,q)\n'.format(k)) 
    for iq,q in enumerate(Q):
        a = S[ik][iq]/S[ik][0]
        with open(ofname, 'a+') as f:
            f.write('{0}\t{1}\n'.format(q, a))
