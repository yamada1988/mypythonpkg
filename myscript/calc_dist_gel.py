#calc_dist.py
import mdtraj as md
import numpy as np
import sys

index_ = int(sys.argv[1])
frame_ = int(sys.argv[2])

Nchain = 50
Nmol = 15841 
d_th = 0.50
fname = 'SYS/DAT/npt_par{0:04d}_nopbc{1:05d}.gro'.format(index_, frame_)
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

ow = df[df['name'] == 'O']
pva_c = df[df['name'] == 'c3']

# inter-chain
outf = 'SYS/DAT/odist_par{0:04d}_chain{1:05d}.dat'.format(index_, frame_)

with open(outf, 'wt') as f:
    f.write('#index\tmol1\tchain2\tatm1\tatm2\tx1\ty1\tz1\tx2\ty2\tz2\tO-O dist\n')

with open(outf, 'a+') as f:
    icount = 0
    for i in range(Nchain+1, Nmol+1):
        ocount = 0
        aow = ow[ow['resSeq']==i]
        pos_aow = np.array(aow[['x','y','z']]) 
        atm_aow = np.array(aow['serial'])
        for j in range(1, Nchain+1):
            if i == j:
                continue
            dn = pva_c[pva_c['resSeq']==j]
            pos_dn = np.array(dn[['x','y','z']])
            atm_dn = np.array(dn['serial'])
            #print(j, pos_dn)
            for ipd, pdn in enumerate(pos_dn):
                dr = np.abs(pos_aow - pdn)
                #print(dr/box, np.round(dr/box))
                dr -= np.round(dr/box)*box
                #print(dr/box)
                d = np.sqrt(np.sum(dr**2))
                #print(d)
                if d <= d_th:
                    icount += 1
                    ocount += 1
                    print('molpair: {0:4d}\t{1:4d}'.format(i,j))
                    print('atmpair: {0:4d}\t{1:4d}\t{2:6.5f}'.format(atm_aow[0], atm_dn[ipd], d))
                    f.write('{0:4d}\t{1:3d}\t{2:4d}\t{3:4d}\t{4:4d}\t{5:7.3f}\t{6:7.3f}\t{7:7.3f}\t{8:7.3f}\t{9:7.3f}\t{10:7.3f}\t{11:6.5f}\n'.format(
                    icount, i, j, atm_aow[0], atm_dn[ipd], pos_aow[0,0], pos_aow[0,1], pos_aow[0,2], pos_dn[ipd,0], pos_dn[ipd,1], pos_dn[ipd,2], d))
        print(i,ocount)
