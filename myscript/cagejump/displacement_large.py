import mdtraj as md
import numpy as np
import sys
import os
import matplotlib.animation as animation
import matplotlib as mpl
mpl.use('Agg')

# Parameters
recdt = 0.10e0 # ps
tint = 10000
dt = recdt * tint #ps
Nchain = 1000
sysname = '../SYS/solution.gro'
xtcname = '../Production/md.xtc'

Nframes = 100000
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


sgn = np.array([0.0 for k in range(3)])
t_max = Nframes
#t_max = 10000
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

nmin = anim_index
nmax = nmin
for num_index in range(nmin, nmax+1):
    #print(num_index)
    positions = []
    it = 0
    for t in md.iterload(xtcname,top=sysname):
        pos = t.xyz
        boxs = t.unitcell_lengths
        t_ = t.time
        maxbox = boxs[0,0]
        for ip,p in enumerate(pos):
            #print(t_[ip])
            if it > Nframes:
                sys.exit()
    
            #if it % 1000000 == 0:
            #    print('t: {0:10.5f}(ns)'.format(it*recdt/1000.0))
            box = boxs[ip,0]
            if maxbox <= box:
                maxbox = box
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
    b = list(zip(*positions))
    
    x = b[0]
    y = b[1]
    z = b[2]
    
    N = len(x)
    print(N)
    #frms = np.array([int(i/tb)*float(tb/float(N)) for i in range(N)])
    frms = np.array([float(i)/float(N) for i in range(N)])
    fname = 'disps/large/disp{0:04d}.dat'.format(num_index)
    with open(fname, 'wt') as f:
        f.write('#time(ns\tdisplacement(nm\n')
        for i in range(1,N-tint):
            x0 = x[i]
            y0 = y[i]
            z0 = z[i]
            x_ = x[i+tint]
            y_ = y[i+tint]
            z_ = z[i+tint]
            dr = np.sqrt(np.sum((x_-x0)**2+(y_-y0)**2+(z_-z0)**2))
            #print(i,dr)
            l = '{0:11.4f}\t{1:7.5f}\n'.format(i*recdt/1000.0,dr)
            f.write(l)
