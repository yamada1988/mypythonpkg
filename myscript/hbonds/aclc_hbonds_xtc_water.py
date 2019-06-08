import mdtraj as md
import numpy as np
import sys
import os
from joblib import Parallel, delayed
from scipy.sparse import lil_matrix, csr_matrix
from numba import jit
import pickle
import cupy as cp
import time

tstart=0
tend=5

def read_info(fname='info.txt'):
    with open(fname, 'rt') as f:
        infodict = {line.split(':')[0]:line.split(':')[1].strip() for line in f}
    return infodict         


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

def make_nlist(pos):
    pos = cp.asarray(pos, dtype=np.float32)    
    t0 = time.time()
    int_list = -1*cp.ones((N, lmax), dtype=np.int32)
    update_list = cp.ElementwiseKernel(
        in_params='int32 data, raw float32 pos, float32 cutoff, float32 L', 
        out_params='raw int32 int_list',
        operation=\
        '''
        int l = 0;
        for (int j = i+1; j < _ind.size(); j++){
             double tmp[3] = {0.0};
             double diff = 0.0;
             for (int k = 0; k<3; k++){ 
               tmp[k] = pos[3*j+k] - pos[3*i+k];
               if (tmp[k] < -0.50*L) {
               tmp[k] += L;
               } else if (tmp[k] > 0.50*L) {
               tmp[k] -= L;
               } 
               diff += powf(tmp[k],2);
             }
             double r = sqrt(diff);
             if ( r > 1.0 *cutoff) {
               continue;
             } else{
                     int ind[] = {i, l};
                     int_list[ind] = j;
                     l += 1;
               }
        }
        ''',
        name='update_list')(data, pos, cutoff, L, int_list)
    
    print('update list:', time.time()-t0)
    return cp.asnumpy(int_list)


def calc_hbond(i, a_pos, d_pos, nlist):
    ls = []
    t0 = time.time()
    pio = a_pos[i-1]
    for j in nlist[i-1]:
        if j == -1:
            break
        pjo = a_pos[j]
        r_joio = pjo - pio
        dr_joio = r_joio - np.round(r_joio/box)*box
        #print(dr/box)
        d = np.sqrt(np.sum(dr_joio**2))
        #print(d)
        if 0.001e0 < d <= d_hbond:
            r_jhio = pio - d_pos[2*j] 
            r_jhjo = pjo - d_pos[2*j] 
            theta = angle_between(r_jhio, r_jhjo) * 180.0/np.pi
            #print((i-1)*Nmon+ipa, (j-1)*Nmon+jpa, d, theta)
            if theta <= theta0 :
                #s[(i-1)*Nmon+ipa,(j-1)*Nmon+jpa] = 1
                ls.append([i-1, 2*j])
            r_jhio = pio - d_pos[2*j+1]
            r_jhio = pio - d_pos[2*j+1]
            theta = angle_between(r_jhio, r_jhjo) * 180.0/np.pi
            if theta <= theta0 :
                ls.append([i-1, 2*j+1])
    t1 = time.time()
    print(it, i, t1-t0)
    if ls != []:
        return ls


# Parameters
Nchain = 138518
Nmon = 1
d_hbond = 0.30 # nm
theta0 = 30 # degree
d_nlist = 2.0 # nm
dt_nlist = 200 # 40ps
lmax = 1500 # max_intnum
sysname = './systemr.gro'
xtcname = './nptr_01000.xtc'
#k = md.load(xtcname, top=sysname)
#Nframes = k.n_frames
Nframes = 10
infodict = read_info()
acname = infodict['acceptor']
dnname1 = infodict['donor1']
dnname2 = infodict['donor2']
if dnname2 == '':
    dnname2 = dnname1
print(acname, dnname1, dnname2)


st = md.load(sysname)
top = st.topology
df, b = top.to_dataframe()
acceptor = df[df['name'] == acname]
donor = df[(df['name'] == dnname1) | (df['name'] == dnname2)]

#print(donor)
a_indexes = np.array(acceptor.index)
d_indexes = np.array(donor.index)
#print(a_indexes)
#print(d_indexes)
N_acc = len(a_indexes)
N_dno = len(d_indexes)
try:
    os.makedirs('./sij_debug')
except :
    pass

sij = []

it = 0
for t in md.iterload(xtcname,top=sysname):
    print('it: {0:5d}'.format(it))
    pos = t.xyz
    box = t.unitcell_lengths[0,0]

    for p in pos:
        t0 = time.time()
        s = lil_matrix((N_acc, N_dno))
        a_pos = p[a_indexes]
        d_pos = p[d_indexes]

        if it >= 10:
            sys.exit()
        if it % dt_nlist == 0:
            nlist = make_nlist(a_pos)

        results = Parallel(n_jobs=-1)([delayed(calc_hbond)(n, a_pos, d_pos, nlist[n]) for n in range(1, Nchain+1)])
        print('##########')
        for rs in results:
            if rs is None:
                continue
            for r in rs:
                s[r[0], r[1]] = 1.0e0
        sij.append(s)
        print(it, s)
        with open('sij_debug/sij_{0:05d}.pickle'.format(it), 'wb') as f:
            pickle.dump(s, f)
        t1 = time.time()
        print(it, t1-t0)
        it += 1
