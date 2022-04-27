import numpy as np
import pickle
import os
import sys

def load_dumps(f):
    obj = []
    while 1:
        try:
            obj.append(pickle.load(f))
        except:
            break
    return obj


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    vector = np.array(vector)
    if np.linalg.norm(vector) <= 0.00010:
        normv = 1.0
    else:
        normv = np.linalg.norm(vector)
    return vector / normv


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


def write_dircinfo(fname, pose):
    outdir = pdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(fname, mode='ab') as f:
        pickle.dump(pose, f)

args = sys.argv
mode = args[1].split('-')[1]
if mode not in ['I', 'i', 'inst', 'Inst', 'ins', 'Ins', 'A', 'a', 'ave', 'Ave']:
    sys.exit('Choose -a or -i.')

pdir = 'npt_r_gpu_pickles/'
posb = pdir + 'box.pickle'
posc = pdir + 'com_pos.pickle'
posd = pdir + 'dirc_pos.pickle'
ilist = pdir + 'intlist.pickle'
with open(posb, 'rb') as f:
    boxes = load_dumps(f)

with open(posc, 'rb') as f:
    com_pos = load_dumps(f)

with open(posd, 'rb') as f:
    dirc_pos = load_dumps(f)

with open(ilist, 'rb') as f:
    intlist = load_dumps(f)


d0 = intlist[0]
del intlist[0]
Frames = len(intlist)
Nchain = len(com_pos[0])
Ncom = len(com_pos[0][0])
#print(d0, Frames)
#print(type(d0))


if mode in ['I', 'i', 'inst', 'Inst', 'ins', 'Ins']:
    outf = pdir+'dirc_instant.pickle'
elif mode in ['A', 'a', 'ave', 'Ave']:
    outf = pdir + 'dirc_averaged.pickle'
with open(outf, mode='wb') as f:
    line = '# d0 = {0:5.4f}'.format(d0)
    pickle.dump(line, f)


for it in range(Frames):
    #print('{0:>5d}/{1:>5d}'.format(it, Frames))
    box = boxes[it]
    order_prm = [[0.0e0 for l in range(Ncom-1) ] for i in range(Nchain)]
    count = [[0.0e0 for l in range(Ncom-1) ] for i in range(Nchain)]
    for il in intlist[it]:
        #print(il)
        i,l,j,h = il
        d_pos = com_pos[it][i][l] - com_pos[it][j][h]
        d_pos -= box*np.trunc(d_pos/(box/2.0))
        d = np.linalg.norm(d_pos)
        #print(i, l, j, h)
        #print(com_pos[it][i][l], com_pos[it][j][h])
        try:
            d1 = dirc_pos[it][i][l]
            d2 = dirc_pos[it][j][h]
        except:
            continue
        t_ = angle_between(d1, d2)
        #print(d, np.cos(t_)**2.0)
        #order_prm[i][l] += (35.0e0 * np.cos(t_)**4.0 - 30.0e0* np.cos(t_)**2.0 +3.0e0)/8.0e0
        #order_prm[j][h] += (35.0e0 * np.cos(t_)**4.0 - 30.0e0* np.cos(t_)**2.0 +3.0e0)/8.0e0
        order_prm[i][l] += 1.0/2.0*(3.0*np.cos(t_)**2.0-1.0)
        order_prm[j][h] += 1.0/2.0*(3.0*np.cos(t_)**2.0-1.0)
        count[i][l] += 1.0
        count[j][h] += 1.0 
        #if (i,l) == (0,0) or (j,h) == (0, 0):
        #    print(i,l,j,h, d, 1.0/2.0*(3.0*np.cos(t_)**2.0-1.0), order_prm[0][0], count[0][0])

    for i in range(Nchain):
        for l in range(Ncom-1):
            #print(i, l, order_prm[i][l], count[i][l])
            try:
                order_prm[i][l] /= count[i][l]
            except ZeroDivisionError:
                order_prm[i][l] = 0.050 
    write_dircinfo(outf, order_prm) 
