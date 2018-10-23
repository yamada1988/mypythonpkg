import sys
import math
import numpy as np
import pickle
from scipy.spatial import distance
from scipy import integrate
import time


#
# This script calculate P(k,q) from compickle.
# Sample:
# python script/calc_Prv.py MD/sample0001_01.pickle 
#

def func(r,v):
    r_box = int(math.trunc(r/dr))
    v_box = int(math.trunc(v/dv))
    if k == 0.0e0 and not q == 0.0e0:
        dummy = P[r_box][v_box]*(2.0e0*r**2)*(v*np.sin(q*v)/(q))*(2.0e0*math.pi)**2
    elif not k == 0.0e0 and q == 0.0e0:
        dummy = P[r_box][v_box]*(r*np.sin(k*r)/(k))*(2.0e0*v**2)*(2.0e0*math.pi)**2
    elif k == 0.0e0 and q == 0.0e0:
        dummy = P[r_box][v_box]*(2.0e0*r**2)*(2.0e0*v**2)*(2.0e0*math.pi)**2
    else:
        dummy = P[r_box][v_box]*(r*np.sin(k*r)/(k))*(v*np.sin(q*v)/(q))*(2.0e0*math.pi)**2
    return dummy 

def phi_u(v, betaM):
    return v**2*math.sqrt(betaM/math.pi)**3*math.exp(-betaM/4.0e0*v**2)

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nens = fname.split('_')[0][-4:]
dirname = fname.split('/')[0]
sysname = fname.split('/')[1].split('_')[0]

pi = math.pi

print('load picklefile:')
t0 = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)
print('time:', time.time() - t0)

tN = len(d['x'])
talpha = 0.010e0
tmax = int(tN*talpha)
tN = tmax
box_size = d['L']
R = d['x']
print(tN, box_size)

# Phisical Parameters
N =  len(R[0])
dt = 0.010e0 #(ps)
L = box_size #(nm) 
rho = float(N/(L**3))
mO = 16.00*1.6611296e-27 #(kg)
mH =  1.008*1.6611296e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)
betaM = 1/(kB*T*Minv)*1.0E+06  # (nm/ps)^-2
print(betaM)

# Space-Velocity Parameters
r_min = 0.0e0
r_max = float(L/1.0e0)
rN = 200
dr = (r_max-r_min)/ float(rN-1)
r_ = np.array([r_min + ir*dr for ir in range(rN)])
vt_min = -0.50e0
vt_max = 1.0e0
vtN = 100
dtv = (vt_max-vt_min)/ float(vtN-1)
vt_ = np.array([vt_min + iv*dtv for iv in range(vtN)])

bnum = 1
tN_b = int(tN/bnum)
totalN = int(N*(N-1)/2)
print('tN_b:',tN_b)

print(r_max)
rm0 = 10.0e0
# Calculate rho(k,q,t)
R = d['x']
V = d['v']
G = [[0.0e0 for iv in range(vtN)] for j in range(rN)]
P = np.zeros((rN, vtN))
g = [0.0e0 for ir in range(rN)]
for it in range(tN_b):
    print('it:', it)
    rvec = np.array(R[it])
    vvec = np.array(V[it])
    # Cannot use distance.pdist due to PBC
    xr = rvec[np.newaxis, :] - rvec[:, np.newaxis]
    xr -= L * np.ceil(xr/L) # unset PBC
    rij = np.sqrt(np.sum(xr**2, axis=2))
    r = rij[np.triu_indices(N, k = 1)]
    rm = np.amin(r)
    if rm < rm0:
        rm0 = rm
    vr = vvec[np.newaxis, :] * vvec[:, np.newaxis]
    vrtheta = np.sum(vr, axis=2)
    vrnorm = np.linalg.norm(vvec, axis=1)
    vr2norm = vrnorm[:,np.newaxis]*vrnorm[np.newaxis,:]
    cos_theta = vrtheta / vr2norm
    p2theta = 3/2*cos_theta**2 - 0.50e0
    t = p2theta[np.triu_indices(N, k = 1)]

    p, rax, tax = np.histogram2d(r, t, bins=(rN, vtN), range=((0.0, r_max), (vt_min, vt_max)))
    P += p

print('sanity check:')
ofname_ = 'DAT/P{0:d}_'.format(int(T))
# normalization
print(np.sum(P), P.shape)
dVr = np.zeros(rN)
dVr = 4.0e0*pi*(r_ +dr/2.0e0)**2*dr
for ir in range(rN):
    for iv in range(vtN):
        G[ir][iv] = P[ir][iv]/(float(N-1)*tN_b)

with open(ofname_+'rt.dat', 'wt') as f:
     for iv in range(vtN):
        for ir in range(rN):
            f.write('{0:6.4f}\t{1:6.4f}\t{2:8.6f}\n'.format(vt_[iv], r_[ir], G[ir][iv]))
        f.write('\n')
fname = 'DAT/S{0:d}_rt.dat'.format(int(T))
Sr = np.zeros(rN)
with open(fname, 'w+') as f:
    for ir in range(rN):
        for iv in range(vtN):
            Sr[ir] += (vt_[iv] + dtv/2.0e0) * G[ir][iv] * dtv
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(r_[ir], Sr[ir]))


with open(ofname_+'r.dat', 'wt') as f:
    g = np.sum(P, axis=1)
    print(rho)
    for ir in range(rN):
        g[ir] /= ((N-1)*tN_b*dVr[ir])
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(r_[ir], g[ir]))

with open(ofname_+'t.dat', 'wt') as f:
    g = np.sum(P, axis=0)
    for iv in range(vtN):
        g[iv] /= (tN_b*N)
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(vt_[iv], g[iv]))


Ginfo = {'r':r_, 'vt':vt_, 'G':G, 'N':N, 'tN':tN_b}
opname = ofname_ + 'rt'.format(T) + '.pickle'
with open(opname, mode='wb') as f:
    pickle.dump(Ginfo, f)
print('dump picklefile time:', time.time() - t0)
