import numpy as np
import sys
import os

prob = []
dth = 0.150
nmin = 1
nmax = 50
dt = 0.00010 #ns
for num_index in range(nmin, nmax+1):
    fname = 'travels_num/water{0:04d}.dat'.format(num_index)
    #fname2 = 'disps/disp{0:04d}.dat'.format(num_index)
    with open(fname, 'rt') as f:
        pairs = [line.split()[2:] for line in f if not line.startswith('#')]
    pairs = [list(map(int, l)) for l in pairs]
    N = len(pairs)
    #print(len(pairs))
    #with open(fname2, 'rt') as f2:
    #    disps = [float(line.split()[1]) for line in f2]   
    #print(len(disps))

    print('num_index:', num_index)
    fname2 = 'break_bond/break_bond{0:04d}.dat'.format(num_index)
    with open(fname2, 'wt') as f2:
        f2.write('#time\tnum\tflag\tpairs\n')

        for it in range(N):
            #print(it)
            if it > 0:
                pairs0 = pairs[it-1]
                try:
                    pairs1 = pairs[it]
                except:
                    break
                hbs = set(pairs0) ^ set(pairs1)
                hbs = list(hbs)
                l = '{0:8.4f}\t{1:d}\t'.format(it*dt, len(hbs))
                if len(hbs) >= 4:
                    l += '1.0\t'
                elif len(hbs) == 3:
                    l += '0.5\t'
                else:
                    l += '0.0\t'
                l += '\t'.join(list(map(str,hbs)))
            else:
                l = '{0:8.4f}\t{1:d}\t'.format(it*dt, 0) 

            l += '\n'
            f2.write(l)
