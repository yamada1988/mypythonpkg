import mdtraj as md
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import matplotlib.cm as cm


args = sys.argv
mode = args[1]

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
            color_[i] = cm.hsv(float(colorid)/N)
    for i in range(N):
         # Make the grid
         try:
             x,y,z= grid[i][:,0], grid[i][:,1], grid[i][:,2]
         except:
             continue
         # Make the direction data for the arrows
         u,v,w = tuple(map(unit_vector, (vec[i][:,0], vec[i][:,1], vec[i][:,2])))
         ax.quiver(x, y, z, u, v, w, length=0.10, normalize=True, cmap = 'Reds')

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

sysgro = './md_run_com.gro'
Nchain = 72
Nmol = 20

t = md.load(sysgro)
pos = t.xyz
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

com = ['' for i in range(Nchain)]
indexes = ['' for i in range(Nchain)]
com_pos = ['' for i in range(Nchain)]
for j in range(Nchain):
    com_ = df[df['name'] == 'c3' ]
    com[j] = com_[com_['resSeq'] == j+1]
    indexes[j] = np.array(com[j].index)
    com_pos[j] = pos[0][indexes[j]]

Ncom = len(com_pos[0])
xyz = ['x', 'y', 'z']

align_pos = ['' for i in range(Nchain)]
for i in range(Nchain):
    align_pos[i] = np.array([[0.0e0 for k in xyz] for j in range(Ncom-1)])
    for l in range(Ncom-1):
        align_pos[i][l] = 0.50e0*(com_pos[i][l+1]-com_pos[i][l])

    # PBC treatment
    for k in range(len(xyz)):
        TorF = com_pos[i][:,k]>box
        #print(xyz[k], TorF)
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            #print(tof_, tof0)
            if tof_ == True:
                com_pos[i][co,k] -= box
            if tof_ != tof0:
                #print(i,co,xyz[k],tof_,tof0)
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
 
        TorF = com_pos[i][:,k]<0.0
        #print(xyz[k], TorF)
        tof0 = TorF[0]
        co = 0
        for tof_ in TorF:
            #print(tof_, tof0)
            if tof_ == True:
                com_pos[i][co,k] += box
            if tof_ != tof0:
                #print(i,co,xyz[k],tof_,tof0)
                align_pos[i][co-1] = np.zeros(3)
                tof0 = tof_
            co += 1
 

if mode in ['Seg', 'seg', 'SEG', 'SEGMENT', 'segment', 's', 'S']:
    Rods = align_pos
    points_pos = ['' for i in range(Nchain)]
    for i in range(Nchain):
        points_pos[i] = np.array([np.array([0.0e0 for k in ['x', 'y', 'z']]) for j in range(Ncom-1)])
        for l in range(Ncom-1):
            points_pos[i][l] = 0.50e0*(com_pos[i][l+1]+com_pos[i][l])
            #print(l, points_pos[i][l], align_pos[i][l])

    d0 = 0.450 #nm
    #ths = [float(it) for it in range(5,90+1, 5)]
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

    #print('crystalline:')
    #for i,c in enumerate(crys_pos):
    #    print(i,len(c))
    
    #print('amorphous:')
    #for i,a in enumerate(amor_pos):
    #    print(i, len(a))

    fig, ax = plot_fig(points_pos, align_pos, box, colormode='mono')
    plt.show()
    sys.exit()
    fig_crys, ax_crys = plot_fig(crys_pos, align_crys, box, fig=None, ax=None, writemode='write', colormode='mono', colorid=0)
    fig_crysamor, ax_crysamor = plot_fig(amor_pos, align_amor, box, fig=fig_crys, ax=ax_crys, writemode='append', colormode='mono', colorid=0.60)
    plt.show()


if mode in ['Rod', 'rod', 'ROD', 'r', 'R', 'ROd', 'rOD', 'rOd']:
    t0 = 30.0
    Rods = [[[]] for i in range(Nchain)]
    for i in range(Nchain):
        k = 0
        ik = 0
        v1 = align_pos[i][k]
        while True:
            if Rods[i][ik] == []:
                Rods[i][ik].append(k)
            v2 = align_pos[i][k+1]
            theta = angle_between(v1, v2)
            #print(i, ik, v1, v2, theta)
            if theta <= t0:
                Rods[i][ik].append(k+1)
                k += 1
            else:
                v1 = align_pos[i][k+1]
                ik += 1
                k += 1
                Rods[i].append([])
            if k+1 == len(align_pos[i]):
                del Rods[i][-1]
                del Rods[i][-1]
                break
    
    for i in range(Nchain):
        for ik in range(len(Rods[i])):
            print(i, ik, Rods[i][ik])
