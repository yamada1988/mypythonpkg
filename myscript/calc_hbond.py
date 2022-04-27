import mdtraj as md
import numpy as np
import sys

index_ = int(sys.argv[1])

fname = 'npt_par{0:04d}_nopbc.gro'.format(index_)
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

acceptor = df[df['name'] == 'oh']
donor = df[df['name'] == 'ho']

print(acceptor)
print(donor)

Nmol = 1611
Nchain = 50
Nstart = 4654
d_hbond = 0.240 # (nm)
# inter-chain
outf = 'hbonds_par{0:04d}_chain.dat'.format(index_)

with open(outf, 'wt') as f:
    f.write('#index\tmol1\tchain2\tatm1\tatm2\tx1\ty1\tz1\tx2\ty2\tz2\tO-H dist\n')

with open(outf, 'a+') as f:
    icount = 0
    for i in range(1, Nchain+1):
        ac = acceptor[acceptor['resSeq']==i]
        pos_ac = np.array(ac[['x','y','z']]) 
        atm_ac = np.array(ac['serial'])
        for j in range(1, Nchain+1):
            if i == j:
                continue
            dn = donor[donor['resSeq']==j]
            pos_dn = np.array(dn[['x','y','z']])
            atm_dn = np.array(dn['serial'])
            #print(j, pos_dn)
            for ipa, pac in enumerate(pos_ac):
                for ipd, pdn in enumerate(pos_dn):
                    dr = np.abs(pac - pdn)
                    #print(dr/box, np.round(dr/box))
                    dr -= np.round(dr/box)*box
                    #print(dr/box)
                    d = np.sqrt(np.sum(dr**2))
                    #print(d)
                    if d <= d_hbond:
                        icount += 1
                        print('molpair: {0:4d}\t{1:4d}'.format(i,j))
                        print('atmpair: {0:4d}\t{1:4d}\t{2:6.5f}'.format(atm_ac[ipa], atm_dn[ipd], d))
                        f.write('{0:4d}\t{1:3d}\t{2:4d}\t{3:4d}\t{4:4d}\t{5:7.3f}\t{6:7.3f}\t{7:7.3f}\t{8:7.3f}\t{9:7.3f}\t{10:7.3f}\t{11:6.5f}\n'.format(
                                icount, i, j, atm_ac[ipa], atm_dn[ipd], pos_ac[ipa,0], pos_ac[ipa,1], pos_ac[ipa,2], pos_dn[ipd,0], pos_dn[ipd,1], pos_dn[ipd,2], d))
