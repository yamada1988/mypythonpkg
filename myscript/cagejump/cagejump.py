import numpy as np
import sys
import os

prob = []
dth = 0.150
nmin = 77
nmax = 1000
N = 100000
dt = 0.0001 #ns
for num_index in range(nmin, nmax+1):
    fname = 'travels_num/water{0:04d}.dat'.format(num_index)
    with open(fname, 'rt') as f:
        pairs = [line.split()[2:] for line in f if not line.startswith('#')]
    pairs = [list(map(int, l)) for l in pairs]
    #print(len(pairs))

    t_interval = 0.0
    flag_jump = 1
    flag_cage = 1
    pairs_cage = ['' for i in range(4)]
    i_cage = 0
    tcages = [['', ''] for i in range(100)]
    for it in range(N):
        #print(it)
        if flag_cage == 1:
            #check state
            if len(pairs[it]) == 4:
                pairs_cage = pairs[it]
                flag_cage = 0
                t0 = it*dt
                #print(tcages,i_cage)
                tcages[i_cage] = [t0,'']
        if flag_cage == 0:
            for pc in pairs_cage:
                if pc in pairs[it]:
                    flag_jump = 1
                    break
                else:
                    flag_jump = 0
            if flag_jump == 0:
                flag_cage = 1
                t1 = it*dt
                tcages[i_cage][1] = t1
                i_cage += 1

    tcages = [l for l in tcages if not '' in l]
    print('num_index:', num_index)
    print('tcage:')

    print(tcages)
    #sys.exit()
    its = []
    for ts in tcages:
        te = ts[1]
        it = int(te/dt)
        its.append(it)

    fname = 'cages/cage{0:04d}.dat'.format(num_index)
    with open(fname, 'wt') as f:
        f.write('#time\thbond\n')
        for it in range(N):
            if it > 0:
                if it in its:
                    l = '{0:7.3f}\t{1:4.2f}'.format(it*dt, 20.0)
                else:
                    l = '{0:7.3f}\t{1:3.2f}'.format(it*dt, 0.0)
            else:
                l = '{0:7.3f}\t{1:3.2f}'.format(it*dt, 0.0)
            l += '\n'
            f.write(l)
