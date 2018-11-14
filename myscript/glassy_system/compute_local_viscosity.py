import h5py
import numpy as np

with h5py.File('MD/local0001_240.h5', 'r') as f:
    lm = f['localMass'].value
    lg = f['localMomentum'].value
    ls = f['localStressTensor'].value
    mg = f['MeshGrid'].value

t = len(lm)
rN = len(lm[0])
dr = mg[0][1]
lv = np.zeros((t,3,rN,rN,rN))
lp = np.zeros((t,rN,rN,rN))
lt = np.zeros((t,3,3,rN,rN,rN))
ldv = np.zeros((t,3,3,rN, rN, rN))
ld = np.zeros((t,rN, rN, rN))
le = np.zeros((t, 3,3,rN, rN, rN))
for it in range(t):
    lv[it] = lg[it]/lm[it]

    for xa in range(3):
        for xb in range(3):
            ls[it][xa][xb] -= lg[it][xa]*lv[it][xb]

    lp[it] = -1.0e0/3.0e0 * np.trace(ls[it])

    for xa in range(3):
        for xb in range(3):
            lt[it][xa][xb] += lp[it]*np.eye(3)[xa][xb]

    ld[it] = np.sum(np.gradient(lv, dr, axis=0))

    for xa in range(3):
        for xb in range(3):
            ldv[it][xa][xb] = 0.50e0*(np.gradient(lv[it][xa], dr, axis=xb)+np.gradient(lv[it][xb], dr,axis=xa)) - 1.0/3.0*ld[it]*np.eye(3)[xa][xb]

    for xa in range(3):
        for xb in range(3):
            le[it][xa][xb] = lt[it][xa][xb] / (2.0*ldv[it][xa][xb])

    print(ls[it,0,0,0,0,0])
