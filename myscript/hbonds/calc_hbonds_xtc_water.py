import mdtraj as md
import numpy as np
import sys
import os
from scipy.sparse import csr_matrix
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
    int_list = -1*cp.ones((Nchain, lmax), dtype=np.int32)
    update_list = cp.ElementwiseKernel(
        in_params='int32 data, raw float32 pos, float32 cutoff, float32 L', 
        out_params='raw int32 int_list',
        operation=\
        '''
        int l = 0;
        for (int j = 1; j < _ind.size(); j++){
             if (i == j){
             continue; 
             }
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
        name='update_list')(data, pos, d_nlist, box, int_list)
    
    print('update list:', time.time()-t0)
    return cp.asnumpy(int_list)


def calc_hbond(a_pos, d_pos, int_list):
    a_pos = cp.asarray(a_pos, dtype=np.float32)    
    d_pos = cp.asarray(d_pos, dtype=np.float32)
    int_list = cp.asarray(int_list, dtype=np.int32)
    h_list = -1*cp.ones((Nchain, 10), dtype=np.int32)
    t0 = time.time()
    calc_hbonds = cp.ElementwiseKernel(
        in_params='int32 data, raw int32 int_list, raw float32 a_pos, raw float32 d_pos, int32 lmax, float32 cutoff, float32 L, float32 rad0', 
        out_params='raw int32 h_list',
        operation=\
        '''
        int m = 0;
        for (int j = 0; j < lmax-1; j++){
             double l = int_list[i*lmax+j];
             double tmp[3] = {0.0};
             double diff = 0.0;
             for (int k = 0; k<3; k++){ 
               tmp[k] = a_pos[3*l+k] - a_pos[3*i+k];
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
                 double tmp1[3] = {0.0};
                 double diff1 = 0.0;
                 double tmp2[3] = {0.0};
                 double diff2 = 0.0;
                 for (int k = 0; k<3; k++){ 
                   tmp1[k] = a_pos[3*i+k] - d_pos[3*2*l+k];
                   if (tmp1[k] < -0.50*L) {
                     tmp1[k] += L;
                   } else if (tmp1[k] > 0.50*L) {  
                     tmp1[k] -= L;
                   }
                 diff1 += powf(tmp1[k],2);
                 }
                 for (int k = 0; k<3; k++){ 
                   tmp2[k] = a_pos[3*l+k] - d_pos[3*2*l+k];
                   if (tmp2[k] < -0.50*L) {
                     tmp2[k] += L;
                   } else if (tmp2[k] > 0.50*L) {  
                     tmp2[k] -= L;
                   }
                 diff2 += powf(tmp2[k],2);
                 }
                 double prod = 0.0;
                 for (int k = 0; k<3; k++){ 
                   prod += tmp1[k]*tmp2[k];
                 }
                 diff1 = sqrt(diff1);
                 diff2 = sqrt(diff2);
                 prod /= (diff1*diff2);
                 if (rad0 <= prod){
                   if (prod <= 1.0){
                     int ind[] = {i, m};
                     h_list[ind] = 2*l;
                     m += 1;
                   }
                 }
                 double tmp3[3] = {0.0};
                 double diff3 = 0.0;
                 double tmp4[3] = {0.0};
                 double diff4 = 0.0;
                 for (int k = 0; k<3; k++){ 
                   tmp3[k] = a_pos[3*i+k] - d_pos[3*2*l+3+k];
                   if (tmp3[k] < -0.50*L) {
                     tmp3[k] += L;
                   } else if (tmp3[k] > 0.50*L) {  
                     tmp3[k] -= L;
                   }
                 diff3 += powf(tmp3[k],2);
                 }
                 for (int k = 0; k<3; k++){ 
                   tmp4[k] = a_pos[3*l+k] - d_pos[3*2*l+3+k];
                   if (tmp4[k] < -0.50*L) {
                     tmp4[k] += L;
                   } else if (tmp4[k] > 0.50*L) {  
                     tmp4[k] -= L;
                   }
                 diff4 += powf(tmp4[k],2);
                 }
                 diff3 = sqrt(diff3);
                 diff4 = sqrt(diff4);
                 prod = 0.0;
                 for (int k = 0; k<3; k++){ 
                   prod += tmp3[k]*tmp4[k];
                 }
                 prod /= (diff3*diff4);
                 if (rad0 <= prod){
                   if (prod <= 1.0){
                     int ind[] = {i, m};
                     h_list[ind] = 2*l+1;
                     m += 1;
                   }
                 }
               }
        }
        ''',
        name='calc_hbonds')(data, int_list, a_pos, d_pos, lmax, d_hbond, box, rad0, h_list)  
    print('calc hbonds:', time.time()-t0)
    return cp.asnumpy(h_list)

def flatten(nested_list):
    return [e for inner_list in nested_list for e in inner_list]


# Parameters
recdt = 0.0040e0 # ps
tint = 25
dt = recdt * tint #ps
Nchain = 138518
Nmon = 1
d_hbond = 0.30 # nm
theta0 = 30 # degree
import math
rad0 = theta0/180.0*math.pi # rad
d_nlist = 1.50 # nm
dt_nlist = 500 # 2.0ps
dt_rec = 500
lmax = 1000 # max_intnum  lmax=1120/2nm
data = cp.arange(Nchain, dtype=np.int32)
sysname = './sys/solution.gro'
xtcname = './MD/npt_r_gpu.xtc'
outdir = 'sij_'+xtcname.split('/')[-1].split('.')[0]
logfile = outdir + '/sij.log'
from datetime import datetime as dat
import datetime
tdatetime = dat.now()
tstr = tdatetime.strftime('%Y-%m-%d-%H:%M:%S')
with open(logfile, 'wt') as f:
    f.write('# sij calculation started at '+ tstr + '\n')
    f.write('# prog%\ttime(ps\ttotbond\tend time\n')
#k = md.load(xtcname, top=sysname)
#Nframes = k.n_frames
Nframes = 100000
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
    os.makedirs(outdir)
except :
    os.remove(outdir+'/sij.pickle')
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

it = 0
tstart = time.time()
tstarttime = dat.now()
for t in md.iterload(xtcname,top=sysname):
    pos = t.xyz
    boxs = t.unitcell_lengths

    for ip,p in enumerate(pos):
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
            nlist = make_nlist(a_pos)

        results = calc_hbond(a_pos, d_pos, nlist)
        print('##########')
        #t0=time.time()
        # 10 times Faster than lil_matrix
        row = np.array([[ir]*len(r[r>=0]) for ir,r in enumerate(results)])
        row = flatten(row)
        col = results[results>=0]
        datas = np.ones(len(col))
        #print(row_data)
        #print(col_data)
        #print([1 for i in range(len(col_data))])
        #for ir,rs in enumerate(results):
        #    for r in rs:
        #        if r >=0:
        #            s[ir, r] = 1.0e0
        s = csr_matrix((datas, (row, col)), shape=(N_acc, N_dno))
        #print('set sij:', time.time()-t0)
        #t0 = time.time()
        with open(outdir+'/sij.pickle', 'ab') as f:
            pickle.dump([recdt*it, s], f, protocol=2)
        #t1 = time.time()
        #print('write picklefile:', t1-t0)

        if it % dt_rec == 0 and it > 0:
            tintertime = dat.now()
            tdelta = tintertime - tstarttime
            per = float(it)/Nframes
            remper = (1.0 - doneper)
            if remper <= 0.0:
                tstr = dat.now().strftime('%Y-%m-%d-%H:%M:%S')
            else:
                tremainsec = tdelta.seconds / (float(it)/Nframes) * remper
                tdelta = datetime.timedelta(seconds=tremainsec)
                tend = tintertime + tdelta
                tstr = tend.strftime('%Y-%m-%d-%H:%M:%S')
            with open(logfile, 'a+') as f:
                f.write('{0:>6.2f}\t{1:>7.2f}\t{2:7d}\t'.format(100.0*per, recdt*it, int(s.sum())) + tstr + '\n')
        
        it += 1
