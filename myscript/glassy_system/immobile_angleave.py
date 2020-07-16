import mdtraj as md
import numpy as np

xtc = 'MD/md/md.xtc'
gro = 'SYS/solution.gro'

t = md.load(xtc,top=gro)
rs = t.xyz
Ls = t.unitcell_lengths
sigma = 0.60 #nm
Nt = 10000
dt = 0.10 #ps
N = t.n_atoms
h_abs = 2.0*np.pi/sigma
rho_i = np.array([0.0E0 for i in range(N)])
for it in range(0,100000,20):
    print("it=",it)
    drs = 0.0E0
    drs_vecsum = np.zeros((N,3))
    for jt in range(Nt):
        L = Ls[it+jt,0]
        drs_vec = rs[it+jt+1] - rs[it+jt]
        drs_vec -= L*np.round(drs_vec/L)
        drs_vecsum += drs_vec
    drs = np.sqrt(np.sum(drs_vecsum**2,axis=1))
    #print(drs)
    rho = 0.0
    for in_atm in range(N): 
        theta = h_abs*drs[in_atm]
        if theta == 0.0:
            rho_i[in_atm] = 1.0E0
        else:
            rho_i[in_atm] = np.sin(theta)/theta
        rho += rho_i[in_atm]/N



    fname = 'DAT/angle_ave/rho_hr{0:05d}.dat'.format(int(it*dt))
    #maxe = np.max(rho_i)
    #mine = np.min(rho_i)
    #de = maxe - mine
    #rho_i -= mine
    de = 1.0E0
    with open(fname, 'w') as f:
        for i in range(N):
            line = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:10.7f}\n".format(rs[it,i,0], rs[it,i,1], rs[it,i,2], rho_i[i]/de)
            f.write(line)
