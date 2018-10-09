import numpy as np
import sys

args = sys.argv
nens_min= int(args[1])
nens_max = int(args[2])
nk_min = int(args[3])
nk_max = int(args[4])
nq_min = int(args[5])
nq_max = int(args[6])
T = int(args[7])

a = [0 for i in range(nens_min, nens_max+1)]
for nk in range(nk_min, nk_max+1):
    for nq in range(nq_min, nq_max+1):
        for ii,i in enumerate(range(nens_min, nens_max+1)):
            ifname = 'DAT/{0:04d}nk{1:02d}nq{2:02d}_{3:03d}.dat'.format(i,nk,nq,T)
            print(ifname)
            with open(ifname, 'rt') as f:
                commentline = [line for line in f if line.startswith('#')]
                f.seek(0)
                a[ii] = np.array([list(map(float, np.array(line.split()))) for line in f if not line.startswith('#')])
        tN = len(a[0])
        ave = [np.array([0.0e0,0.0e0,0.0e0]) for it in range(tN)] 
        for ii,i in enumerate(range(nens_min, nens_max+1)):
            ave += a[ii]
        ave /= float(nens_max - nens_min + 1)
        ofname = 'DAT/Cnk{0:02d}nq{1:02d}_{2:03d}.dat'.format(nk, nq, T)
        with open(ofname, 'wt') as f:
            for line in commentline:
                f.write(line)
        with open(ofname, 'a+') as f:
            for it in range(tN):
                f.write('{0:9.5f}\t{1:6.5f}\t{2:6.5f}\n'.format(ave[it][0], ave[it][1], ave[it][2]))
