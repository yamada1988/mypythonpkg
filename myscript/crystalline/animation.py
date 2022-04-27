import mdtraj as md
import numpy as np
import sys
import os
import pickle
import time
import matplotlib.animation as animation

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
    int_list = np.array(int_list)
    dist_pos = a1 - a2
    dist_pos -= box * np.trunc(dist_pos/(box/2.0))
    d = np.sqrt(np.sum(dist_pos**2, axis=1))
    d_index = (0.001e0 < d ) & (d < d0)
    #print(d_index)
    #print(int_list)
    #print(int_list[d_index])
    ids = [id0 for id0,d0 in enumerate(d_index) if d0 == True]
    #print(ids)
    for ip in range(2):
        did = 2*int_list[d_index] + ip
        #print(did)
        #print(a1, d_pos0[did])
        r_jhio = a1 - d_pos0[did]
        r_jhio -= box * np.trunc(r_jhio/(box/2.0))
        r_jhjo = a2[ids] - d_pos0[did]
        r_jhjo -= box * np.trunc(r_jhjo/(box/2.0))
        #print(r_jhio, r_jhjo)
        r_jhio = np.array([unit_vector(r) for r in r_jhio])
        r_jhjo = np.array([unit_vector(r) for r in r_jhjo])
        #if it == 16:
        #    print(i, d, d_index, r_jhio, r_jhjo)
        try:
            thetas = np.sum(r_jhio*r_jhjo, axis=1)
        except ValueError:
            break 
        #print(thetas)
        theta1s = np.arccos(np.clip(thetas,-1.0,1.0)) * 180.0/np.pi
        #print(theta1s)
        t_index = (theta1s <= theta0 ) | (180.0-theta0 <= theta1s)
        #print(t_index)
        dtid = int_list[d_index][t_index]
        #print(dtid)
        for id_ in dtid:
            latm.append(2*id_+ip)
            lij.append(id_)

    #print('calc hbonds:', time.time()-t0)
    return latm, lij

def flatten(nested_list):
    return [e for inner_list in nested_list for e in inner_list]


# Parameters
recdt = 0.10e0 # ps
tint = 100
dt = recdt * tint #ps
Nchain = 1000
d_hbond = 0.350 # nm
theta0 = 30 # degree
import math
rad0 = theta0/180.0*math.pi # rad
sysname = '../SYS/solution.gro'
xtcname = '../Production/md.xtc'
outdir = 'sij_'+xtcname.split('/')[-1].split('.')[0] + '_skip/'
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
Nframes = 1000000
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


it = 0
tstart = time.time()
tstarttime = dat.now()
sgn = np.array([0.0 for k in range(3)])
t_max = 50000

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm



for num_index in range(1, 1+1):
    positions = []
    for t in md.iterload(xtcname,top=sysname):
        pos = t.xyz
        boxs = t.unitcell_lengths
        t_ = t.time
        maxbox = boxs[0,0]
        for ip,p in enumerate(pos):
            #print(t_[ip])
            if it > Nframes:
                sys.exit()
            if it % tint != 0:
                it += 1
                continue
    
            print('it: {0:5d}'.format(it))
            box = boxs[ip,0]
            if maxbox <= box:
                maxbox = box
            t0 = time.time()
            a_pos = p[a_indexes]
            a = a_pos[num_index-1]
            if positions != []:
                dr =  a__ - a 
                for k in range(3):
                    if (dr[k]>= box*0.50):
                        sgn[k] += 1.0
                    elif (dr[k]<= -box*0.50):
                        sgn[k] -= 1.0
                #print(a,dr,sgn)
            a_ = a+sgn*box
            a__ = a
            positions.append(a_)
            #print(a,a_,sgn)
            
            it += 1
        if it >= t_max:
            break
    
    positions = np.array(positions)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Set label
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    # Set Box
    box_ = maxbox
    
    b = list(zip(*positions))
    x = b[0]
    y = b[1]
    z = b[2]
    
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    zmin = np.min(z)
    zmax = np.max(z)
    
    ax.set_xlim(xmin-0.10, xmax+0.10)
    ax.set_ylim(ymin-0.10, ymax+0.10)
    ax.set_zlim(zmin-0.10, zmax+0.10)

tb = 25

N = len(x)
print(N)
def animate(i):
    x_ = x[i]
    y_ = y[i]
    z_ = z[i]
    n = int(i/tb)*float(tb/N)
    print(i,n)
    cm_ = cm.hsv(n)
    ax.scatter(x_, y_, z_,s=40,alpha=0.8, marker="o",c=cm_)
    if i%tb == 0:
        ax.set_title("time:{0:6.2f} (ps)".format(i*dt), 
             loc="left",
             fontdict={"fontsize":   18,
                       "fontweight": "bold"})

# 表示
ani =animation.FuncAnimation(fig, animate,frames=N,repeat=False) 
plt.show()
plt.close()
