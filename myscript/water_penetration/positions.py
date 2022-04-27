import numpy as np
import sys
import os

args = sys.argv
phi = float(args[1])

L = 32.83307  
rho = 700.0

mp = 44.0
nump = 200
Np = 500
Mp = (mp+2.0)*nump*Np

Mw = phi/(1.0-phi)*Mp
Vtarget = 0.50*Mw/(6.23*10**23)*10**-3/rho
dL = Vtarget/(L*L*10**(-18))/(10**-9) #nm
Nw = int(Mw/18.0)
nw1 = int(Nw/2.0)
nw2 = Nw - nw1
print(dL, nw1)

Lmax1 = 20.0
Lmin1 = Lmax1 - dL

Lmin2 = 52.883
Lmax2 = Lmin2 + dL

fname = 'positions.dat'
try:
    os.remove(fname)
except:
    pass

for i in range(1, nw1+1):
    x = L * np.random.rand()  
    y = L * np.random.rand()
    z = (Lmax1 - Lmin1) * np.random.rand() + Lmin1 
    #print(x,y,z)
    with open(fname, 'a+') as f:
        f.write('{0:8.3f}\t{1:8.3f}\t{2:8.3f}\n'.format(x,y,z))


for i in range(nw1+1, Nw):
    x = L * np.random.rand()  
    y = L * np.random.rand()
    z = (Lmax2 - Lmin2) * np.random.rand() + Lmin2
    with open(fname, 'a+') as f:
        f.write('{0:8.3f}\t{1:8.3f}\t{2:8.3f}\n'.format(x,y,z))
