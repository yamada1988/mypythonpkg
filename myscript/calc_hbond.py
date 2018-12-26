import mdtraj as md
import numpy as np

fname = 'npt_par_low05_0001.gro'
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
#print(donor)

Nchain = 50

# inter-chain
outf = 'hbonds_inter_chain.dat'
with open(outf, 'wt') as f:
    f.write('#chain1\tchain2\tatm1\tatm2\tO-O dist\n')

with open(outf, 'a+') as f:
    for i in range(1, Nchain+1):
        ac = acceptor[acceptor['resSeq']==i]
        pos_ac = np.array(ac[['x','y','z']]) 
        atm_ac = np.array(ac['serial'])
        for j in range(i+1, Nchain+1):
            dn = acceptor[acceptor['resSeq']==j]
            pos_dn = np.array(dn[['x','y','z']])
            #print(j, pos_dn)
            for ipa, pac in enumerate(pos_ac):
                for ipd, pdn in enumerate(pos_dn):
                    dr = (pac - pdn)
#                    print(dr/box, np.trunc(dr/box))
                    dr -= np.trunc(dr/box)*box
                    d = np.sqrt(np.sum(dr**2))
                    if d <= 0.30:
                        print('molpair: {0:4d}\t{1:4d}'.format(i,j))
                        print('atmpair: {0:4d}\t{1:4d}'.format(atm_ac[ipa], atm_ac[ipd]), d)
                        f.write('{0:2d}\t{1:2d}\t{2:4d}\t{3:4d}\t{4:6.5f}\n'.format(i, j, atm_ac[ipa], atm_ac[ipd], d))
