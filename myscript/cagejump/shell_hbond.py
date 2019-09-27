import numpy as np
import sys
import os

shells = []
nmin = anim_index
nmax = nmin
dt = 0.00010 #ns
for num_index in range(nmin, nmax+1):
    fname = 'travels_num/water{0:04d}.dat'.format(num_index)
    with open(fname, 'rt') as f:
        pairs = [line.split()[2:] for line in f if not line.startswith('#')]
    pairs = [list(map(int, l)) for l in pairs]
    N = len(pairs)

    print('num_index:', num_index)
    fname2 = 'shell_bond/shell_bond{0:04d}.dat'.format(num_index)
    with open(fname2, 'wt') as f2:
        f2.write('#time\tnum\tshell\tpairs\n')

    shell_index = 0
    index_pairs = []
    for it in range(N):
        if it >= 0:
            pairs1 = pairs[it]
            if len(pairs1) == 4:
                if pairs1 not in shells:
                    shell_index += 1
                    index_pairs.append([shell_index, pairs1])
                shells.append(pairs1)

                for ip in index_pairs:
                    if pairs1 == ip[1]: 
                        ind = ip[0]
            if len(index_pairs) == 0: 
                l = '{0:8.4f}\t{1:>3d}\t'.format(it*dt, 0)
            else:
                l = '{0:8.4f}\t{1:>3d}\t'.format(it*dt, ind)
                l += '\t'.join(list(map(str,pairs1)))

            l += '\n'
            with open(fname2, 'at') as f2:
                f2.write(l)
