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
    return v**2*math.sqrt(betaM/math.pi)**3*np.exp(-betaM*v**2)

args = sys.argv
fname = args[1]
T = float(fname.split('_')[1].split('.')[0])
nens = fname.split('_')[0][-4:]
dirname = fname.split('/')[0]
sysname = fname.split('/')[1].split('_')[0]

pi = math.pi

print('load picklefile:')
t = time.time()
with open(fname, mode='rb') as f:
    d = pickle.load(f)
print('time:', time.time() - t)

tN = len(d['x'])
talpha = 1.0e0
tmax = int(tN*talpha)
tN = tmax
box_size = d['L']
R = d['x']
print(tN, box_size)

# Phisical Parameters
N =  len(R[0])
dt = 1.0e0 #(ps)
L = box_size #(nm) 
rho = float(N/(L**3))
mO = 16.00e-27 #(kg)
mH =  1.008e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)
betaM = 1/(kB*T*Minv)*1.0E+06  # (nm/ps)^-2
print(betaM)

# Space-Velocity Parameters
r_min = 0.0e0
r_max = float(L/2.0e0)
rN = 200
dr = (r_max-r_min)/ float(rN)
r_ = np.array([r_min + ir*dr for ir in range(rN)])
v_min = 0.0e0
v_max = 4.0e0
vN = 500
dv = (v_max-v_min)/ float(vN)
v_ = np.array([v_min + iv*dv for iv in range(vN)])

bnum = 1
tN_b = int(tN/bnum)
totalN = int(N*(N-1)/2)
print('tN_b:',tN_b)

print(r_max)
rm0 = 10.0e0
# Calculate rho(k,q,t)
R = d['x']
V = d['v']
G = [[0.0e0 for iv in range(vN)] for j in range(rN)]
P = np.zeros((rN, vN))
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
    vr = vvec[np.newaxis, :] - vvec[:, np.newaxis]
    vij = np.sqrt(np.sum(vr**2, axis=2))
    v = vij[np.triu_indices(N, k = 1)]

    p, rax, vax = np.histogram2d(r, v, bins=(rN, vN), range=((0.0, r_max), (0.0, v_max)))
    P += p

print('sanity check:')
ofname_ = 'DAT/P{0:d}_'.format(int(T))
# normalization
print(np.sum(P), P.shape)
dVr = np.zeros(rN)
dVv = np.zeros(vN)
dVr = 4.0e0*pi*(r_ +dr/2.0e0)**2*dr
dVv = 4.0e0*pi*(v_ +dv/2.0e0)**2*dv
for ir in range(rN):
    for iv in range(vN):
        G[ir][iv] = P[ir][iv]/(float(N-1)*tN_b*dVr[ir]*2.0e0)
with open(ofname_+'rv.dat', 'wt') as f:
     for iv in range(vN):
        for ir in range(rN):
            f.write('{0:6.4f}\t{1:6.4f}\t{2:8.6f}\n'.format(v_[iv], r_[ir], G[ir][iv]))
        f.write('\n')

with open(ofname_+'r.dat', 'wt') as f:
    g = np.sum(P, axis=1)
    print(rho)
    for ir in range(rN):
        g[ir] /= ((N-1)*tN_b*dVr[ir]*2.0e0)
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(r_[ir], g[ir]))

with open(ofname_+'v.dat', 'wt') as f:
    g = np.sum(P, axis=0)
    for iv in range(vN):
        g[iv] /= (tN_b*N)
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(v_[iv], g[iv]))

Ginfo = {'r':r_, 'v':v_, 'G':G, 'N':N, 'tN':tN_b}
opname = ofname_ + 'rv'.format(T) + '.pickle'
with open(opname, mode='wb') as f:
    pickle.dump(Ginfo, f)
print('dump picklefile time:', time.time() - t)
