import mdtraj as md
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import matplotlib.cm as cm


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
    return 180.0e0/np.pi * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def judge(x, threshold):
    if x <= threshold or x >= 180.0-threshold:
        return 1
    else:
        return 0


# Plot
def plot_fig(grid, vec, box_, fig=None, ax=None, writemode='write', colormode='grad', colorid=0):
    N = len(grid)
    for i in range(N):
        grid[i] = np.array(grid[i])
        vec[i] = np.array(vec[i])
    if writemode == 'write':
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Set label
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_zlabel("Z-axis")
        # Set Box
        ax.set_xlim(0.0, box_)
        ax.set_ylim(0.0, box_)
        ax.set_zlim(0.0, box_) 
    
    # Set colormode
    color_ = ['' for i in range(N)]
    if colormode == 'grad':
        for i in range(N):
            color_[i] = cm.hsv(float(i)/N)
    elif colormode == 'mono':
        for i in range(N):
            color_[i] = cm.hsv(float(colorid))

    for i in range(N):
         # Make the grid
         try:
             x,y,z= grid[i][:,0], grid[i][:,1], grid[i][:,2]
         except:
             continue
         # Make the direction data for the arrows
         u,v,w = tuple(map(unit_vector, (vec[i][:,0], vec[i][:,1], vec[i][:,2])))
         ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True, color=color_[i])

    return fig, ax


def remove_empty(data):
    while True:
        try:
            data.remove([])
        except:
            break
    return data


def make_gridvec(ls, pos, align):
    g = [[] for i in range(Nchain)]
    v = [[] for i in range(Nchain)]
    for l in ls:
        #print(l, pos[l[0]][l[1]])
        g[l[0]].append(pos[l[0]][l[1]])
        v[l[0]].append(align[l[0]][l[1]])

    return g, v

sysgro = './md_run_director.gro'
Nchain = 72
comgro = sysgro.split('_director.gro')[0] + '_com.gro'
d0 = 0.450 #nm

with open(sysgro, 'rt') as f:
    align_pos = np.array([np.array(list(map(float, line.split()[6:9]))) for line in f if len(line) >= 48])
    f.seek(0)
    lines = f.readlines()

Ncom = 0
for l in lines:
    if ' 1PE' in l:
        Ncom += 1

align_pos = align_pos.reshape(Nchain, Ncom, 3)

with open(comgro, 'rt') as f:
    com_pos = np.array([np.array(list(map(float, line.split()[3:6]))) for line in f if len(line) >= 28])
    f.seek(0)
    boxes = f.readlines()[-1].split()

points_pos = ['' for i in range(Nchain)]
com_pos = com_pos.reshape(Nchain, Ncom+1, 3)
for i in range(Nchain):
    points_pos[i] = ['' for k in range(Ncom)]
    for k in range(Ncom):
       points_pos[i][k] = com_pos[i][k]

points_pos = np.array(points_pos)
box = float(boxes[0])

xyz = ['x', 'y', 'z']

# PBC treatment for align_pos
for i in range(Nchain):
    for k in range(len(xyz)):
        TorF = com_pos[i][:,k]>box
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            if tof_ == True:
                points_pos[i][co,k] -= box
            if tof_ != tof0:
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
            if co == Ncom:
                break
 
        TorF = com_pos[i][:,k]<0.0
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            if tof_ == True:
                points_pos[i][co,k] += box
            if tof_ != tof0:
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
            if co == Ncom:
                break
 


crys_list = []
crys_pos = [[] for i in range(Nchain)]
align_crys = [[] for i in range(Nchain)]
amor_pos = [[] for i in range(Nchain)]
align_amor = [[] for i in range(Nchain)]
ths = [20.0]
success = [0]
count = 0
for i in range(Nchain):
    #ls = [s for s in range(Ncom-1)]
    for j in range(i+1,Nchain):
        for l0 in range(len(points_pos[i])):
            for l in range(len(points_pos[i])):
                d_pos = points_pos[i][l0] - points_pos[j][l]
                d = np.linalg.norm(d_pos)
                if d < d0:
                    #print(i,l0, j,l)
                    #if l0 in ls:
                    #    ls.remove(l0)
                    #print(i,ls)
                    t_ = angle_between(align_pos[i][l0], align_pos[j][l])
                    #print(i,l0, j,l, d, align_pos[i][l0], align_pos[j][l], t_)
                    #print(i,l0,j,l,judge(t_, th), d,t_)
                    if np.linalg.norm(align_pos[i][l0]) > 0.0 and np.linalg.norm(align_pos[i][l0]) > 0.0:
                        count += 1
                        for ith,th in enumerate(ths):
                            if np.linalg.norm(align_pos[i][l0]) > 0.0 and np.linalg.norm(align_pos[i][l0]) > 0.0:  
                                success[ith] += judge(t_,th)
                                if t_ <= th or th > 180.0-th:
                                    if [i, l0] not in crys_list:
                                        crys_list.append([i,l0])
                                    if [j,l] not in crys_list:
                                        crys_list.append([j,l])


amor_list = [[i,l] for i in range(Nchain) for l in range(Ncom-1)]
for list_ in crys_list:
    try:
        amor_list.remove(list_)
    except:
        continue


crys_pos, align_crys = make_gridvec(crys_list, points_pos, align_pos)
amor_pos, align_amor = make_gridvec(amor_list, points_pos, align_pos)


for ith, th in enumerate(ths):
   print(th, float(success[ith])/float(count), count)


fig, ax = plot_fig(points_pos, align_pos, box, colormode='grad')
plt.show()
sys.exit()
#fig_crys, ax_crys = plot_fig(crys_pos, align_crys, box, fig=None, ax=None, writemode='write', colormode='mono', colorid=0.0)
#fig_amor, ax_amor = plot_fig(amor_pos, align_amor, box, fig=None, ax=None, writemode='write', colormode='mono', colorid=0.50)
#fig_crysamor, ax_crysamor = plot_fig(amor_pos, align_amor, box, fig=fig_crys, ax=ax_crys, writemode='append', colormode='mono', colorid=0.60)
#plt.show()
