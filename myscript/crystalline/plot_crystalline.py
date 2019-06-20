import numpy as np
import pickle
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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


def write_dircinfo(fname, pose):
    outdir = pdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(fname, mode='ab') as f:
        pickle.dump(pose, f)

args = sys.argv
fmode = args[1].split('-')[1]
if fmode not in ['I', 'i', 'inst', 'Inst', 'ins', 'Ins', 'A', 'a', 'ave', 'Ave']:
    sys.exit('Choose -a or -i.')

pdir = 'npt_r_gpu_pickles/'
posb = pdir + 'box.pickle'
posc = pdir + 'com_pos.pickle'
posd = pdir + 'dirc_pos.pickle'

if fmode in ['I', 'i', 'inst', 'Inst', 'ins', 'Ins']:
    posi = pdir + 'dirc_instant.pickle'
    with open(posi, 'rb') as f:
        order_prm = load_dumps(f)
elif fmode in ['A', 'a', 'ave', 'Ave']:
    posi = pdir + 'dirc_averaged.pickle'
    with open(posi, 'rb') as f:
        order_prm = load_dumps(f)

with open(posb, 'rb') as f:
    boxes = load_dumps(f)
with open(posc, 'rb') as f:
    com_pos = load_dumps(f)
with open(posd, 'rb') as f:
    dirc_pos = load_dumps(f)

#print(order_prm[0])
del order_prm[0]
Frames = len(order_prm)
Nchain = len(com_pos[0])
Ncom = len(com_pos[0][0])

# Plot
def plot_fig(grid, vec, o_prm, box_, cmode='ordered'):
    Nc = len(grid)
    Nm = len(grid[0])
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

    for i in range(Nc):
        for l in range(Nm-1):
            #print(i,l, o_prm[i][l])
            grid_ = np.array(grid[i][l])
            vec_ = np.array(vec[i][l]) 
            # Set colormode
            if cmode == 'ordered':
                c_ = cm.spectral(1.0-o_prm[i][l])
        
            # Make the grid
            try:
                x,y,z= grid_[0], grid_[1], grid_[2]
            except:
                continue
            # Make the direction data for the arrows
            u,v,w = vec_[0], vec_[1], vec_[2]
            ax.quiver(x, y, z, u, v, w, length=0.2, normalize=True, color=c_)

    return fig, ax

for it in range(Frames):
    for i in range(Nchain):
        for l in range(Ncom-1):
            if order_prm[it][i][l] < 0.10:
                order_prm[it][i][l] = 0.10e0

for it in range(Frames):
    fig, ax = plot_fig(com_pos[it], dirc_pos[it], order_prm[it], boxes[it])
    figname = pdir + 'order_parameter{0:04d}.eps'.format(it)
    plt.savefig(figname)
    plt.show()
