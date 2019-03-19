import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import mdtraj as md
import numpy as np
import sys

colordict = matplotlib.colors.cnames
colors =  list(colordict.keys())
print(colors)

# Parameters
Nchain = 50
d_hbond = 0.240 # (nm)

fname = 'pbc_1000.gro'
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]


#print(df)
tmp = df[df['resName']=='PVA20']
backbones = tmp[tmp['name']=='c3'].reset_index()
acceptor = df[df['name'] == 'oh']
donor = df[df['name'] == 'ho']

print(backbones)
#print(acceptor)
#print(donor)

sns.set_palette("hls", Nchain)
cf = 'clusters.dat'

# Hydrogen bonding between inter-chain
with open(cf, 'rt') as f:
    pos_clust = [[float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]

bf = 'cluster_hbonds.dat'
with open(bf, 'rt') as f:
    c_atms = [[int(line.split()[0]), int(line.split()[4]), int(line.split()[5])] for line in f if not line.startswith('#')]
print(pos_clust)

# Plot backbones
fig = pyplot.figure(figsize=(12,8))
ax = Axes3D(fig)

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-0.0, 1.0*box)
ax.set_ylim(-0.0, 1.0*box)
ax.set_zlim(-0.0, 1.0*box)

# Plot cluster atoms
Nc = len(pos_clust)
for i in range(Nc):
    for ca in c_atms:
        if ca[0] == i+1:
            x1 = backbones[(backbones['serial']==ca[1]-2)]['x'].values
            y1 = backbones[(backbones['serial']==ca[1]-2)]['y'].values
            z1 = backbones[(backbones['serial']==ca[1]-2)]['z'].values
            x2 = backbones[(backbones['serial']==ca[2]-3)]['x'].values
            y2 = backbones[(backbones['serial']==ca[2]-3)]['y'].values
            z2 = backbones[(backbones['serial']==ca[2]-3)]['z'].values
            ax.plot(x1,y1,z1, "o", color=colors[i%50], ms=4.0,markeredgewidth=1, markeredgecolor='black')
            ax.plot(x2,y2,z2, "o", color=colors[i%50], ms=4.0,markeredgewidth=1, markeredgecolor='black')

for ca in c_atms:
    backbones.drop(backbones.index[backbones.serial==ca[1]-2], inplace=True)
    backbones.drop(backbones.index[backbones.serial==ca[2]-3], inplace=True)


xs = ['' for i in range(Nchain)]
ys = ['' for i in range(Nchain)]
zs = ['' for i in range(Nchain)]
for i in range(1, Nchain+1):
    print(backbones[(backbones['resSeq']==i)])
    xs[i-1] = backbones[(backbones['resSeq']==i)]['x'].values
    ys[i-1] = backbones[(backbones['resSeq']==i)]['y'].values
    zs[i-1] = backbones[(backbones['resSeq']==i)]['z'].values
    ax.plot(xs[i-1], ys[i-1], zs[i-1], ".", color='black', markersize=2.0, markeredgewidth=1)


# Plot clusters
for i in range(Nc):
    x = pos_clust[i][0]
    y = pos_clust[i][1]
    z = pos_clust[i][2]
    print('index:',i, x, y, z)
    #ax.plot([x], [y], [z], "o", color=colors[i%50], ms=25, alpha=0.50)
    #ax.text(x, y, z, '{0:2d}'.format(i+1), color='white', fontsize=16)


ofname = 'npt.eps'
#pyplot.savefig(ofname)

pyplot.show()
