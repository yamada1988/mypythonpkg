import mdtraj as md
import numpy as np
import sys
import os
from joblib import Parallel, delayed
from scipy.sparse import lil_matrix, csr_matrix
import pickle
import time

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

def calc_hbond(i, a_pos, d_pos):
    ls = []
    for j in range(1, Nchain+1):
        t0 = time.time()
        if j == i:
            #continue
            pass
        for ipa, pio in enumerate(a_pos[(i-1)*Nmon:i*Nmon]):
            for jpa, pjo in enumerate(a_pos[(j-1)*Nmon:j*Nmon]):
                if ipa == jpa:
                    continue
                r_joio = pjo - pio
                dr_joio = r_joio - np.round(r_joio/box)*box
                #print(dr/box)
                d = np.sqrt(np.sum(dr_joio**2))
                #print(d)
                if 0.001e0 < d <= d_hbond:
                    r_jhio = pio - d_pos[(j-1)*Nmon+jpa] 
                    r_jhjo = pjo - d_pos[(j-1)*Nmon+jpa] 
                    theta = angle_between(r_jhio, r_jhjo) * 180.0/np.pi
                    #print((i-1)*Nmon+ipa, (j-1)*Nmon+jpa, d, theta)
                    if theta <= theta0 :
                        #s[(i-1)*Nmon+ipa,(j-1)*Nmon+jpa] = 1
                        ls.append([(i-1)*Nmon+ipa, (j-1)*Nmon+jpa])
        t1 = time.time()
        #print(it, i, j, t1-t0)
    if ls != []:
        return ls


# Parameters
Nchain = 100
Nmon = 100
d_hbond = 0.30 # (nm)
theta0 = 30 # (degree)
sysname = './systemr.gro'
xtcname = './npt_r_gpu.xtc'
k = md.load(xtcname, top=sysname)
Nframes = k.n_frames

st = md.load(sysname)
top = st.topology
df, b = top.to_dataframe()
acceptor = df[df['name'] == 'ho']
donor = df[df['name'] == 'oh']
a_indexes = np.array(acceptor.index)
d_indexes = np.array(donor.index)
#print(a_indexes)
#print(d_indexes)
reses = np.array([i+1 for i in range(Nmon) for k in range(100)])
try:
    os.makedirs('./sij')
except :
    pass

sij = []

it = 0
for t in md.iterload(xtcname,top=sysname):
    print('it: {0:5d}'.format(it))
    pos = t.xyz
    boxs = t.unitcell_lengths

    #print(pos.shape)
    for ip,p in enumerate(pos):
        box = boxs[ip,0]
        if it < tstart or it >= tend:
            print(it)
            it += 1
            continue
        t0 = time.time()
        s = lil_matrix((Nchain*Nmon, Nchain*Nmon))
        a_pos = p[a_indexes]
        d_pos = p[d_indexes]
        results = Parallel(n_jobs=-1)([delayed(calc_hbond)(n, a_pos, d_pos) for n in range(1, Nchain+1)])
        print('##########')
        for rs in results:
            if rs is None:
                continue
            for r in rs:
                s[r[0], r[1]] = 1.0e0
        sij.append(s)
        print(it, s)
        with open('sij/sij_{0:05d}.pickle'.format(it), 'wb') as f:
            pickle.dump(s, f)
        t1 = time.time()
        print(it, t1-t0)
        it += 1
