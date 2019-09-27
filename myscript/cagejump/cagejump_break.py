import numpy as np
import sys
import os

prob = []
dth = 0.150
nth = 4
nmin = 101
nmax = 700
N = 50000
dt = 0.0010 #ns
for num_index in range(nmin, nmax+1):
    fname = 'break_bond/break_bond{0:04d}.dat'.format(num_index)
    fname2 = 'disps/disp{0:04d}.dat'.format(num_index)
    print(num_index)
    with open(fname, 'rt') as f:
        break_num = [int(line.split()[1]) for line in f if not line.startswith('#')]
    #print(len(pairs))
    with open(fname2, 'rt') as f2:
        disps = [float(line.split()[1]) for line in f2 if not line.startswith('#')]   
    #print(len(disps))

    its = []
    for it in range(N):
        if break_num[it] >= nth:
            #print(it)
            its.append(int(it*1.0)+1)
            d_max = np.max(np.array(disps[int(it/1.0)-1:int(it/1.0)+2]))
            #print(te, disps[it-1:it+2])
            #print(te, d_max)
            if d_max > dth:
                prob.append(1.0)
            else:
                prob.append(0.0)
    fname = 'cages_break/cage_break{0:04d}.dat'.format(num_index)
    with open(fname, 'wt') as f:
        f.write('#time\tdispl\thbond\n')
        for it in range(int(N/1.0)):
            if it > 0:
                l = '{0:7.3f}\t{1:5.4f}'.format(it*dt*1.0, disps[it-1])
            else:
                l = '{0:7.3f}\t{1:5.4f}'.format(it*dt*1.0, 0.0)
            if it in its: 
                l += '\t   {0:3.2f}'.format(dth)
            #if disps[it-1] > dth:
            #    l += '\t   {0:3.2f}'.format(dth+0.010)
            l += '\n'
            f.write(l)

print(sum(prob)/len(prob))
