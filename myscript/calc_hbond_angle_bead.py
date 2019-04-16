import mdtraj as md
import numpy as np
import sys

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



fname = '../systemlllll_0050.gro'
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

acceptor = df[df['name'] == 'ho']
donor = df[df['name'] == 'oh']

print(acceptor)
print(donor)

Nchain = 1
d_hbond = 0.320 # (nm)
theta0 = 30 # (degree)
# inter-chain
outf = 'hbonds_chain.dat'

with open(outf, 'wt') as f:
    f.write('#index\tmol1\tchain2\tatm1\tatm2\tx1\ty1\tz1\tx2\ty2\tz2\tr_OO\t angle\n')

with open(outf, 'a+') as f:
    icount = 0
    for i in range(1, Nchain+1):
        ac1 = acceptor[acceptor['resSeq']==i]
        pos_o1 = np.array(ac1[['x','y','z']]) 
        atm_o1 = np.array(ac1.index)
        dn1 = donor[donor['resSeq']==i]
        pos_h1 = np.array(dn1[['x','y','z']])
        atm_h1 = np.array(dn1.index)
        for j in range(1, Nchain+1):
            #if i == j:
            #    continue
            acj = acceptor[acceptor['resSeq']==j]
            pos_o2 = np.array(acj[['x','y','z']])
            atm_o2 = np.array(acj.index)
            dn2 = donor[donor['resSeq']==j]
            pos_h2 = np.array(dn2[['x','y','z']])
            atm_h2 = np.array(dn2.index)
            #print(j, pos_dn)
            for ipa, po1 in enumerate(pos_o1):
                for jpa, po2 in enumerate(pos_o2):
                    r_o2o1 = po2 - po1
                    #print(dr/box, np.round(dr/box))
                    dr_o2o1 -= np.round(r_o2o1/box)*box
                    #print(dr/box)
                    d = np.sqrt(np.sum(dr_o2o1**2))
                    #print(d)
                    if d <= d_hbond:
                        r_h1o1 = pos_h1[ipa] - pos_o1[ipa] 
                        theta = angle_between(r_h1o1, r_o2o1) * 180.0/np.pi
                        if theta <= theta0:
                            icount += 1
                            print('molpair: {0:4d}\t{1:4d}'.format(i,j))
                            print('atmpair: {0:4d}\t{1:4d}\t{2:6.5f}\t{3:3.2f}'.format(atm_h1[ipa], atm_o2[jpa], d, theta))
                            print(pos_o2[jpa], pos_h1[ipa])
                            print(pos_o1[ipa]) 
                            f.write('{0:4d}\t{1:3d}\t{2:3d}\t{3:6d}\t{4:6d}\t{5:7.3f}\t{6:7.3f}\t{7:7.3f}\t{8:7.3f}\t{9:7.3f}\t{10:7.3f}\t{11:6.5f}\t{12:3.2f}\n'.format(
                                    icount, i, j, atm_h1[ipa]+1, atm_o2[jpa]+1, pos_h1[ipa][0], pos_h2[ipa][1], pos_h1[ipa][2], pos_o2[jpa][0], pos_o2[jpa][1], pos_o2[jpa][2], d, theta))
