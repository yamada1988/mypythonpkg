import sys
import math
import numpy as np
import pickle
from scipy.spatial import distance
from scipy import integrate
import time


#
# This script calculate P(x,y,z,ngp) from compickle.
# Sample:
# python script/calc_Prv.py MD/sample0001_01.pickle 
#

def f(ingpu, ix, iy, iz, zp):
    return ngp_[ingp]*P[ix][iy][iz][ingp] / zp

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
mO = 16.00*1.6611296e-27 #(kg)
mH =  1.008*1.6611296e-27 #(kg)
Minv = 1/(mO+mH+mH) #(kg^-1)
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
ngp_max = 5.0e0
ngp_min = -5.0e0
ngpN = 1000
dngp = (ngp_max - ngp_min) / float(ngpN)
ngp_ = np.array([ngp_min + ingp*dngp for ingp in range(ngpN)])
ngp_s = ngp_ + dngp*0.50e0

# Calculate P(x,y,z,u,v,w)
R = np.array(d['x'])
V = d['v']
# unset PBC 
R -= np.trunc(R/L)*L
R += np.round((0.50e0*L-R)/L)*L

P = np.zeros((rN, rN, rN, ngpN))
g = [0.0e0 for ir in range(rN)]
for it in range(tN):
    print('it:', it)
    rvec = np.array(R[it])
    vvec = np.array(V[it]) 
    ngps = (vvec**2 - 1.0e0/betaM)
    ngps = np.sum(ngps, axis=1)
    print(ngps)
    rngp = np.c_[rvec,ngps]
    p, rnax = np.histogramdd(rngp, bins=(rN, rN, rN, ngpN), range=((0.0, r_max), (0.0, r_max), (0.0, r_max), (ngp_min, ngp_max)), normed=True)
    P += p

# normalization
P /= (tN)
rnax = np.array(rnax)
print(rnax)
print(rnax.shape)
print(P.shape)
print('sanity check:')
print('Psum:', np.sum(P*dr*dr*dr*dngp))

ofname_ = 'DAT/Prngp{0:d}.dat'.format(int(T))
with open(ofname_, 'wt') as f:
    for ix in range(rN):
        for iy in range(rN):
            for iz in range(rN):
                for ingp in range(ngpN):
                    f.write('{0:6.4f}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:8.7f}\n'.format(r_[ix], r_[iy], r_[iz], ngp_[ingp], P[ix][iy][iz][ingp]))
                f.write('\n')

for ix in range(rN):        
    ofname_ = 'DAT/manyfiles/Ngpr{0:02d}_{1:d}.dat'.format(ix, int(T))
    with open(ofname_, 'wt') as f:
        for iy in range(rN):
            for iz in range(rN):
                np = 0.0e0
                zp = 0.0e0
                for ingp in range(ngpN):
                    zp += P[ix][iy][iz][ingp]*dngp
                fp = ngp_s * P[ix][iy][iz]/zp * dngp
                np = integrate.simps(fp, range(ngpN))
                f.write('{0:6.4f}\t{1:6.4f}\t{2:8.7f}\n'.format(r_[iy], r_[iz], np))
            f.write('\n')

sys.exit()
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
