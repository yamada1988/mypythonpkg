import mdtraj as md
import numpy as np
import sys
import os
from scipy.sparse import csr_matrix
from numba import jit
import pickle
from joblib import Parallel, delayed
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


def make_nlist(pos, box, d0):
    t0 = time.time()
    pos = np.array(pos, dtype=np.float32)
    dist_pos = pos[:,np.newaxis] - pos[np.newaxis,:]
    dist_pos -= box * np.trunc(dist_pos/(box/2.0))
    #print(pos.shape, dist_pos.shape)
    d = np.sqrt(np.sum(dist_pos**2, axis=2))
    torfs = (0.0 < d) & (d < d0)
    int_list = [list(np.where(t)[0]) for t in torfs]
    print('update list:', time.time()-t0)
    #print(int_list)
    return int_list


def calc_hbond(i, a1, a2, d_pos0, d0, int_list, box):
    latm = []
    lij = []
    a1 = np.array([a1])
    dist_pos = a1 - a2
    dist_pos -= box * np.trunc(dist_pos/(box/2.0))
    d = np.sqrt(np.sum(dist_pos**2, axis=1))
    for jd,d_ in enumerate(d):
        #print(jd,d_)
        if  0.001e0 < d_ <= d0:
            ja = int_list[jd]
            r_jhio = a1 - d_pos0[2*ja]
            r_jhio -= box * np.trunc(r_jhio/(box/2.0))
            r_jhjo = a2[jd] - d_pos0[2*ja]
            r_jhjo -= box * np.trunc(r_jhjo/(box/2.0))
            theta1 = angle_between(r_jhio, r_jhjo) * 180.0/np.pi
            if theta1 <= theta0 or 180.0-theta0 <= theta1:
                latm.append(2*ja)
                lij.append(ja)

            r_jhio = a1 - d_pos0[2*ja+1]
            r_jhio -= box * np.trunc(r_jhio/(box/2.0))
            r_jhjo = a2[jd] - d_pos0[2*ja+1]
            r_jhjo -= box * np.trunc(r_jhjo/(box/2.0))
            theta2 = angle_between(r_jhio, r_jhjo) * 180.0/np.pi
            if theta2 <= theta0 or 180.0-theta0 <= theta2:
                latm.append(2*ja+1)
                if ja not in lij:
                    lij.append(ja)

    #print('calc hbonds:', time.time()-t0)
    return latm, lij

def flatten(nested_list):
    return [e for inner_list in nested_list for e in inner_list]


# Parameters
recdt = 0.10e0 # ps
tint = 1
dt = recdt * tint #ps
Nchain = 5000
d_hbond = 0.320 # nm
theta0 = 30 # degree
import math
rad0 = theta0/180.0*math.pi # rad
d_nlist = 1.20 # nm
dt_nlist = 10 # 5.0ps
dt_rec = 50
sysname = '../SYS/solution.gro'
xtcname = '../MD/md.xtc'
outdir = 'sij_'+xtcname.split('/')[-1].split('.')[0]
logfile = outdir + '/sij.log'
try:
    os.makedirs('./'+outdir)
except :
    pass

from datetime import datetime as dat
import datetime
tdatetime = dat.now()
tstr = tdatetime.strftime('%Y-%m-%d-%H:%M:%S')
with open(logfile, 'wt') as f:
    f.write('# sij calculation started at '+ tstr + '\n')
    f.write('# prog%\ttime(ps\ttotbond\tend time\n')
#k = md.load(xtcname, top=sysname)
#Nframes = k.n_frames
Nframes = 50000
acname = 'OW1'
dnname1 = 'H2'
dnname2 = 'HW3'
print(acname, dnname1, dnname2)


st = md.load(sysname)
top = st.topology
df, b = top.to_dataframe()
acceptor = df[df['name'] == acname]
donor = df[(df['name'] == dnname1) | (df['name'] == dnname2)]

a_indexes = np.array(acceptor.index)
d_indexes = np.array(donor.index)
#print(d_indexes)
N_acc = len(a_indexes)
N_dno = len(d_indexes)
try:
    os.makedirs(outdir)
except :
    pass
try:
    os.remove(outdir+'/sij.pickle')
except:
    pass
try:
    os.remove(outdir+'/satm.pickle')
except:
    pass

with open(outdir+'/sij.pickle', 'wb') as f:
    line = '''
This picklefile contains time series of hydrogen-bond informations for each sites, sij, as csr_matrix.
Import scipy.sparse module if you want to read.
Example for load sij.pickle:
====================
def load_dumps(f):
    obj = []
    while 1:
        try:
            obj.append(pickle.load(f))
        except:
            break
    return obj
if __name__ == '__main__':
    with open('sij.pickle', 'rb') as f:
        data = load_dumps(f)
    del data[0]
====================
'''
    pickle.dump(line, f)

with open(outdir+'/satm.pickle', 'wb') as f:
    pickle.dump(line, f)




it = 0
tstart = time.time()
tstarttime = dat.now()

for t in md.iterload(xtcname,top=sysname):
    pos = t.xyz
    boxs = t.unitcell_lengths
    t_ = t.time
    for ip,p in enumerate(pos):
        print(t_[ip])
        if it > Nframes:
            sys.exit()
        if it % tint != 0:
            it += 1
            continue

        print('it: {0:5d}'.format(it))
        box = boxs[ip,0]
        t0 = time.time()
        a_pos = p[a_indexes]
        d_pos = p[d_indexes]

        if it % dt_nlist == 0:
            nlist = make_nlist(a_pos, box, d_nlist)

        t0 = time.time()
        results = [calc_hbond(i, a_pos[i], a_pos[il], d_pos, d_hbond, il, box) for i,il in enumerate(nlist)]
        results = list(zip(*results))
        results_atm = results[0]
        results_ij = results[1]
        print(time.time()-t0)
        #for ir,r in enumerate(results):
        #    print(ir,r)
      
        print('##########')
        #t0=time.time()
        # 10 times Faster than lil_matrix
        row = np.array([[ir]*len(r) for ir,r in enumerate(results_atm)])
        row = flatten(row)
        col = flatten(results_atm)
        datas = np.ones(len(col))
        #print(row_data)
        #print(col_data)
        #print([1 for i in range(len(col_data))])
        #for ir,rs in enumerate(results):
        #    for r in rs:
        #        if r >=0:
        #            s[ir, r] = 1.0e0
        s_atm = csr_matrix((datas, (row, col)), shape=(N_acc, N_dno), dtype=np.int32)
        with open(outdir+'/satm.pickle', 'ab') as f:
            pickle.dump([t_[ip], s_atm], f, protocol=2)

        row = np.array([[ir]*len(r) for ir,r in enumerate(results_ij)])
        row = flatten(row)
        col = flatten(results_ij)
        datas = np.ones(len(col))
        s_ij = csr_matrix((datas, (row, col)), shape=(Nchain, Nchain), dtype=np.int32)
        with open(outdir+'/sij.pickle', 'ab') as f:
            pickle.dump([t_[ip], s_ij], f, protocol=2)


        if it % dt_rec == 0 and it > 0:
            tintertime = dat.now()
            tdelta = tintertime - tstarttime
            per = float(it)/Nframes
            remper = (1.0 - per)
            if remper <= 0.0:
                tstr = dat.now().strftime('%Y-%m-%d-%H:%M:%S')
            else:
                tremainsec = tdelta.seconds / (float(it)/Nframes) * remper
                tdelta = datetime.timedelta(seconds=tremainsec)
                tend = tintertime + tdelta
                tstr = tend.strftime('%Y-%m-%d-%H:%M:%S')
            with open(logfile, 'a+') as f:
                f.write('{0:>6.2f}\t{1:>7.2f}\t{2:7d}\t'.format(100.0*per, recdt*it, int(s_atm.sum())) + tstr + '\n')

        
        it += 1
