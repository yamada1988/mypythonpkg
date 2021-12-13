import numpy as np
import sys

args = sys.argv
index = args[1]

temperature = 483.0 #K
kBT = 2.479*temperature/298.0 #kJ/mol
beta = 1.0/kBT

fname1 = 'DAT/energy{0:04d}_S0.dat'.format(int(index))
fname2 = 'DAT/energy{0:04d}_S0B0.dat'.format(int(index))

es_S0 = np.loadtxt(fname1)
es_S0B0 = np.loadtxt(fname2)

es_diff = es_S0 - es_S0B0

fname3 = 'DAT/energy{0:04d}_B0.dat'.format(int(index))

es_B0 = np.loadtxt(fname3)

es_exact = es_B0 + es_diff

#print(es_exact)

ofname1 = 'DAT/energy{0:04d}_exact.dat'.format(int(index))

with open(ofname1, 'wt') as f:
    f.write('# step\tenergy(kJ/mol)\n')
    N = len(es_exact)
    for i in range(N):
        line = '{0:04d}\t{1:12.6f}\n'.format(i, float(es_exact[i][1]))
        f.write(line)


fname4 = 'DAT/energy{0:04d}_ghost.dat'.format(int(index))

es_ghost = np.loadtxt(fname4)

es_diff = es_ghost - es_exact

ofname2 = 'DAT/energy{0:04d}_diff.dat'.format(int(index))

with open(ofname2, 'wt') as f:
   f.write('# step\t (W-U)(kJ/mol)\n')
   N = len(es_exact)
   for i in range(N):
       line = '{0:04d}\t{1:12.6f}\n'.format(i, float(es_diff[i][1]))
       f.write(line)

expbDU = np.exp(beta*es_diff[:,1])
ZexpbDU = np.sum(expbDU)

expbDU /= ZexpbDU

ofname3 = 'DAT/SysWght{0:04d}'.format(int(index))
with open(ofname3, 'wt') as f:
    f.write('# step\t exp(beta(W-U)/<exp(beta(W-U)>)\n')
    N = len(expbDU)
    for i in range(N):
        line = '{0:04d}\t{1:e}\n'.format(i, expbDU[i])
        f.write(line)
