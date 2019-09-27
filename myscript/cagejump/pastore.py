import mdtraj as md
import numpy as np
import sys
import os

# Parameters
recdt = 0.10e0 # ps
t_delta = 30
dt = recdt #ps
Nchain = 1000
sysname = '../SYS/solution.gro'
xtcname = '../Production/md.xtc'

print('t_delta:', t_delta)

Nframes = 100000
N = Nframes
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
    
    pos = np.array(positions)
    positions = np.array(list(zip(*pos)))
    print(np.shape(positions)) 

    #print(np.shape(positions[0:40]))   
    frms = np.array([float(i)/float(N) for i in range(N)])
    fname = 'pastores/pastore{0:04d}.dat'.format(num_index)
    with open(fname, 'wt') as f:
        f.write('#time(ns\tflctuation\n')
    for i in range(N):
        if i > t_delta and i < N-t_delta-1:
            #print(positions[i-t_delta:i+t_delta+1])
            x = positions[0][i-t_delta:i+t_delta+1]
            y = positions[1][i-t_delta:i+t_delta+1]
            z = positions[2][i-t_delta:i+t_delta+1]
            x_ave = np.mean(x)
            y_ave = np.mean(y)
            z_ave = np.mean(z)
            x_diff = np.var(x)
            y_diff = np.var(y)
            z_diff = np.var(z)
            #print(x)
            #print(i,x_diff)
            dr = np.sum(x_diff+y_diff+z_diff)
            #print(i,dr)
            l = '{0:11.4f}\t{1:7.5f}\n'.format(i*recdt/1000.0,dr)
        else:
            l = '{0:11.4f}\t{1:7.5f}\n'.format(i*recdt/1000.0, 0.0)
        with open(fname, 'a+') as f:
            f.write(l)
