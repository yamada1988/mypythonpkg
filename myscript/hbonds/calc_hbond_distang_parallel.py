import mdtraj as md
import numpy as np
import sys
from joblib import Parallel, delayed

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
    t1 = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) * 180.0/np.pi
    #v3_u = np.cross(v1_u, v2_u)
    #if v3_u[2] < 0.0:
    #    t1 = 360.0E0 - t1
    return t1



fname = './MD/calc_hbond.gro'
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

# PVA, donor_o = acceptor
acceptor = df[df['name'] == 'oh']
donor_h = df[df['name'] == 'ho']

print(acceptor)
print(donor_h)

Nchain = 500
d_hbond = 0.320 # (nm)
theta0 = 30 # (degree)
# inter-chain
outf = 'hbonds_chain.dat'

with open(outf, 'wt') as f:
    f.write('#index\tmol1\tmol2\tatm1\tatm2\tx1\ty1\tz1\tx2\ty2\tz2\tr_OO\t angle\n')

def calc_hbond(i):
    ls = []
    aci = acceptor[acceptor['resSeq']==i]
    pos_io = np.array(aci[['x','y','z']]) 
    atm_io = np.array(aci.index)
    #dni = donor[donor['resSeq']==i]
    #pos_ih = np.array(dni[['x','y','z']])
    #atm_ih = np.array(dni.index)
    for j in range(1, Nchain+1):
        dnj = acceptor[acceptor['resSeq']==j]
        pos_jo = np.array(dnj[['x','y','z']])
        atm_jo = np.array(dnj.index)
        dhj = donor_h[donor_h['resSeq']==j]
        pos_jh = np.array(dhj[['x','y','z']])
        atm_jh = np.array(dhj.index)
        #print(j, pos_dn)
        for ipa, pio in enumerate(pos_io):
            for jpa, pjo in enumerate(pos_jo):
                if ipa == jpa:
                    continue
                r_DA = pjo - pio
                #print(d= r_joio - (r_/box, np.round(dr/box))
                dr_DA = r_DA - np.round(r_DA/box)*box
                #print('dr:', dr_DA)
                d = np.sqrt(np.sum(dr_DA**2))
                #print(d)
                theta = 0
                if d <= d_hbond:
                    r_HA = pos_jh[jpa] - pos_io[ipa]
                    r_DH = pos_jo[jpa] - pos_jh[jpa]
                    theta = angle_between(r_HA, r_DH)
                    if theta <= theta0 :
                        l = '{0:3d}\t{1:3d}\t{2:6d}\t{3:6d}\t{4:7.3f}\t{5:7.3f}\t{6:7.3f}\t{7:7.3f}\t{8:7.3f}\t{9:7.3f}\t{10:6.5f}\t{11:3.2f}\n'.format(
                            i, j, atm_io[ipa]+1, atm_jh[jpa]+1, pos_io[ipa][0], pos_io[ipa][1], pos_io[ipa][2], pos_jh[jpa][0], pos_jh[jpa][1], pos_jh[jpa][2], d, theta)
                        ls.append(l)
                    #if ipa < 20:
                    #    print(ipa, jpa, atm_io[ipa]+1, atm_jh[jpa]+1, pos_io[ipa], pos_jh[jpa], pos_jo[jpa])
                    #    print('dist:', d, 'angle:', theta)
    return ls


results = Parallel(n_jobs=-1)([delayed(calc_hbond)(n) for n in range(1, Nchain+1)])
icount = 0
for ls in results:
    for l in ls:
        icount += 1
        line = '{0:4d}\t'.format(icount) + l
        with open(outf, 'a+') as f:
            f.write(line)
