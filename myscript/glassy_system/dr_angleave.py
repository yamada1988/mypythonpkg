import mdtraj as md
import numpy as np

xtc = 'MD/md/md.xtc'
gro = 'SYS/solution.gro'

t = md.load(xtc,top=gro)
rs = t.xyz
Ls = t.unitcell_lengths
Nt = 50000
dt = 0.10 #ps
N = t.n_atoms
rho_i = np.array([0.0E0 for i in range(N)])
for it in range(0,100000,100):
    print("it=",it)
    drs = 0.0E0
    drs_vecsum = np.zeros((N,3))
    for jt in range(Nt):
        L = Ls[it+jt,0]
        drs_vec = rs[it+jt+1] - rs[it+jt]
        drs_vec -= L*np.round(drs_vec/L)
        drs_vecsum += drs_vec
    drs = np.sqrt(np.sum(drs_vecsum**2,axis=1))


    fname = 'DAT/drs/dr{0:05d}.dat'.format(int(it*dt))
    with open(fname, 'w') as f:
        for i in range(N):
            line = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:10.7f}\n".format(rs[it,i,0], rs[it,i,1], rs[it,i,2], drs[i])
            f.write(line)
