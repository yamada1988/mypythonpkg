import numpy as np
import sys

target = 100
dr = 0.040 #nm
Lstd = 12.0
gr = np.array([0.0e0 for i in range(int(Lstd/2.0e0*1/dr))])
for i in range(1,100):
    fname = 'SYS/com{0:04d}.gro'.format(i)
    with open(fname, 'rt') as f:
        Rs = [[float(line.split()[-3]), float(line.split()[-2]), float(line.split()[-1])] for line in f if len(line) >= 40]
        f.seek(0)
        L = float(f.readlines()[-1].split()[0])
    print(Rs)
    print(L)
    rho = 10000/(L**3)
    Lhalf = L/2.0e0
    Nmax = int(Lhalf/dr)
    r0 = Rs[target-1]
    r0 = np.array(r0)
    rslv = Rs[target:]
    npos = 10000 - 100
    for rs in rslv:
        rs = np.array(rs)
        r = rs - r0
        #print(r)
        #print(np.round(r/L)) 
        r -= np.round(r/L)*L
        #print(r)
        rnorm = np.linalg.norm(r)
        index = int(rnorm/dr)
        #print(rnorm, index)
        if index <= Nmax:
            gr[index] += 1
    print(i)
for ir in range(Nmax):
    dv = 4.0e0/3.0e0 * np.pi * (ir**3-(ir-1)**3)*dr**3
    gr[ir] /= rho * dv * 100.0e0

with open('rdfs.xvg', 'wt') as f:
    f.write('# r \t gr\n')
with open('rdfs.xvg', 'a+') as f:
    for ir in range(Nmax):
        f.write('{0:6.3f}\t{1:6.4f}\n'.format(ir*dr, gr[ir]))
