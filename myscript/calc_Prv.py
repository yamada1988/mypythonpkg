import sys
import math
import numpy as np
import pickle
from scipy.spatial import distance
from scipy import integrate
import time


#
# This script calculate P(x,y,z,u,v,w) from compickle.
# Sample:
# python script/calc_Prv.py MD/sample0001_01.pickle 
#

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
talpha = 0.00010e0
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
mO = 16.00*1.6611296e-27 #(kg)
mH =  1.008*1.6611296e-27 #(kg)
Minv = 1/(mO + mH + mH) #(kg^-1)
kB = 1.3801e-23 #(J K^-1)
betaM = 1/(kB*T*Minv)*1.0E+06  # (nm/ps)^-2
print(betaM)

# Space-Velocity Parameters
r_min = 0.0e0
r_max = L
rN = 16
dr = (r_max-r_min)/ float(rN)
r_ = np.array([r_min + ir*dr for ir in range(rN)])
v_0 = math.sqrt(betaM)
v_max = 4.0e0
v_min = -4.0e0
vN = 16
dv = (v_max - v_min) / float(vN)
v_ = np.array([v_min + iv*dv for iv in range(vN)])
vs = v_ + dv*0.50e0

# Calculate P(x,y,z,u,v,w)
R = np.array(d['x'])
V = d['v']
# unset PBC 
R -= np.trunc(R/L)*L
R += np.round((0.50e0*L-R)/L)*L

G = [[0.0e0 for iv in range(vN)] for j in range(rN)]
P = np.zeros((rN, rN, rN, vN, vN, vN))
g = [0.0e0 for ir in range(rN)]
for it in range(tN):
    print('it:', it)
    rvec = np.array(R[it])
    vvec = np.array(V[it]) 
    rv = np.hstack((rvec,vvec))
    p, rvax = np.histogramdd(rv, bins=(rN, rN, rN, vN, vN, vN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max), (v_min, v_max), (v_min, v_max), (v_min, v_max)))
    P += p

# normalization
P /= (tN*N)
rvax = np.array(rvax)
print(rvax)
print(rvax.shape)
print(P.shape)
print('sanity check:')
print('ZP:', np.sum(P))

#sys.exit()
with open(ofname_+'rv.dat', 'wt') as f:
     for iv in range(vN):
        for ir in range(rN):
            f.write('{0:6.4f}\t{1:6.4f}\t{2:8.6f}\n'.format(v_[iv], r_[ir], G[ir][iv]))
        f.write('\n')

for ir in range(rN):
    fname = 'DAT/manyfiles/P{0:d}_r{1:03d}.dat'.format(int(T), ir)
    with open(fname, 'wt') as f:
        for iv in range(vN):
            f.write('{0:6.4f}\t{1:8.6f}\t{2:8.6f}\t{3:8.6f}\n'.format(v_[iv], G[ir][iv], phi(vs[iv], betaM), G[ir][iv]/phi(vs[iv], betaM)))


with open(ofname_+'r.dat', 'wt') as f:
    g = np.sum(P, axis=1)
    print(rho)
    for ir in range(rN):
        g[ir] /= (0.50e0*N*(N-1)*tN_b)
        f.write('{0:6.4f}\t{1:8.6f}\n'.format(r_[ir], g[ir]))

with open(ofname_+'v.dat', 'wt') as f:
    g = np.sum(P*dr, axis=0)/ np.sum(P*dr*dv)
    for iv in range(vN):
        f.write('{0:6.4f}\t{1:8.6f}\t{2:8.6f}\t{3:8.6f}\n'.format(v_[iv], g[iv], phi(vs[iv], betaM), g[iv]/phi(vs[iv], betaM)))

Ginfo = {'r':r_, 'v':v_, 'G':G, 'N':N, 'tN':tN_b}
opname = ofname_ + 'rv'.format(T) + '.pickle'
with open(opname, mode='wb') as f:
    pickle.dump(Ginfo, f)
print('dump picklefile time:', time.time() - t)
