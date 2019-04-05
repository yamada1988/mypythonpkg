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



fname = 'npt.gro'
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

Nchain = 100
d_hbond = 0.240 # (nm)
theta0 = 30 # (degree)
# inter-chain
outf = 'hbonds_chain.dat'

with open(outf, 'wt') as f:
    f.write('#index\tmol1\tchain2\tatm1\tatm2\tx1\ty1\tz1\tx2\ty2\tz2\tO-H dist\n')

with open(outf, 'a+') as f:
    icount = 0
    for i in range(1, Nchain+1):
        aci = acceptor[acceptor['resSeq']==i]
        pos_io = np.array(aci[['x','y','z']]) 
        atm_io = np.array(aci['serial'])
        dni = donor[donor['resSeq']==i]
        pos_ih = np.array(dni[['x','y','z']])
        atm_ih = np.array(dni['serial'])
        for j in range(1, Nchain+1):
            if i == j:
                continue
            acj = acceptor[acceptor['resSeq']==j]
            pos_jo = np.array(acj[['x','y','z']])
            atm_jo = np.array(acj['serial'])
            dnj = donor[donor['resSeq']==j]
            pos_jh = np.array(dnj[['x','y','z']])
            atm_jh = np.array(dnj['serial'])
            #print(j, pos_dn)
            for ipa, pio in enumerate(pos_io):
                for jpa, pjo in enumerate(pos_jo):
                    r_joio = pjo - pio
                    dr_joio = np.abs(pio - pjo)
                    #print(dr/box, np.round(dr/box))
                    dr_joio -= np.round(dr_joio/box)*box
                    #print(dr/box)
                    d = np.sqrt(np.sum(dr_joio**2))
                    #print(d)
                    if d <= d_hbond:
                        r_ihio = ipa - pos_ih[ipa] 
                        theta = angle_between(r_ihio, r_joio) * 180.0/np.pi
                        if theta <= theta0:
                            icount += 1
                            print('molpair: {0:4d}\t{1:4d}'.format(i,j))
                            print('atmpair: {0:4d}\t{1:4d}\t{2:6.5f}\t{3:3.2f}'.format(atm_ih[ipa], atm_jo[jpa], d, theta))
                            f.write('{0:4d}\t{1:3d}\t{2:4d}\t{3:4d}\t{4:4d}\t{5:7.3f}\t{6:7.3f}\t{7:7.3f}\t{8:7.3f}\t{9:7.3f}\t{10:7.3f}\t{11:6.5f}\n'.format(
                                    icount, i, j, atm_ih[ipa], atm_jo[jpa], pos_ih[ipa][0], pos_ih[ipa][1], pos_ih[ipa][2], pos_jo[jpa][0], pos_jo[jpa][1], pos_jo[jpa][2], d))
