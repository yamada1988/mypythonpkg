##############################
# Script Calculating  Stress Tensor with Block Averaging Method

# Reference:
# UNDERSTANDING MOLECULAR SIMULATION From Algorithms to Applications
# Daan Frenkel & Berend Smit, ACADEMIC PRESS, 1996



##############################
import numpy as np
import sys
from math import sqrt

args = sys.argv
fname = args[1]

V = 16.784**3.0E0 #(nm^3)
kBT = 2.479 #(kJ/mol)
alpha = kBT/V #(kJ/mol/nm^3)
alpha = alpha * 1.602*10**(6.0E0) #(Pa)


Ps = np.loadtxt(fname, dtype='float', comments=['#', '@'])
Ps = Ps.T
t = Ps[0]
P = Ps[1:]


print(t)

dt = t[1]-t[0]
N = len(P[0])
print(N)
index0 = [0, 10**2, 10**4, 10**6, 10**8]
t0 = [i*dt for i in index0]
bins =   [dt,10*dt,10**2*dt, 10**4*dt, 10**6*dt]
ni = len(index0)+1
ibins = [1, 10, 10**2, 10**4, 10**6]
index = index0 + [N]

Ps = [['' for k in range(1, ni)] for l in range(3) ]
ts = ['' for k in range(1, ni)]
ts[0] = t
nptxt = ['' for k in range(1, ni)]
for k in range(1, ni):
    for l in range(3):
        print('#####\nl,k={0:d}, {1:d}\n######'.format(l, k))
        b = ibins[k]
        Ps[l][k] = np.array([sum(P[l][i*b:(i+1)*b])/b for i in range(int(N/b))])
        ts[k] = [i*bins[k] for i in range(int(N/b))]
    nptxt[k] = [ts[k], Ps[0][k], Ps[1][k], Ps[2][k]]
    nptxt[k] = np.array(list(zip(*nptxt[k])))
    fname = 'sigma_{0:02d}.xvg'.format(k)
    np.savetxt(fname, nptxt[k], delimiter='\t', header='Pxy \tPxz \tPyz (bar)', fmt='%.10e')
