import numpy as np
import sys
import os

hist = []
dth = 0.110
nmin =  1
nmax = 500
N = 50000 - 1
dt = 0.0010 #ns
times = [[] for it in range(N)]
for num_index in range(nmin, nmax+1):
    fname = 'disps/disp{0:04d}.dat'.format(num_index)
    try:
        with open(fname, 'rt') as f:
            drs = [float(line.split()[1]) for line in f if not line.startswith('#')]
        #print(fname, N,len(drs))
    except:
        continue
    for it in range(N):
        if drs[it] >= dth:
            times[it].append(num_index)

fname = 'disps/tjump_rec.dat'
with open(fname, 'wt') as f:
    f.write('#time(ns\tindex\n')
    for it in range(N):
        #print(it*dt, times[it])
        l = '{0:10.4f}\t'.format(it*dt)
        if len(times[it]) != 0:
            l += '\t'.join(list(map(str,times[it])))
        l += '\n'
        f.write(l)
