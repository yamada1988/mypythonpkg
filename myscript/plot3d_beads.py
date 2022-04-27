import matplotlib
matplotlib.use('Qt4Agg')
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
Nchain = 100
d_hbond = 0.240 # (nm)

ifname = 'pbc_r5000.gro'
index = ifname.split('_')[1].split('.')[0]
t = md.load(ifname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]
Natm = 140500 

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]


#print(df)
tmp = df[df['resName']=='PVA20']
backbones = tmp[tmp['name']=='c3'].reset_index()
chain = np.array(backbones['index'].values)
chain += 1
chain = chain.tolist() 
print(chain)
acceptor = df[df['name'] == 'oh']
donor = df[df['name'] == 'ho']

print(backbones)
#print(acceptor)
#print(donor)

sns.set_palette("hls", Nchain)
cf = 'clusters_' +index+ '.dat'

# Hydrogen bonding between inter-chain
with open(cf, 'rt') as f:
    pos_clust = [[float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]

bf = 'cluster_hbonds_' + index+ '.dat'
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
            x1 = backbones[(backbones['index']==ca[1]-3)]['x'].values
            y1 = backbones[(backbones['index']==ca[1]-3)]['y'].values
            z1 = backbones[(backbones['index']==ca[1]-3)]['z'].values
            x2 = backbones[(backbones['index']==ca[2]-4)]['x'].values
            y2 = backbones[(backbones['index']==ca[2]-4)]['y'].values
            z2 = backbones[(backbones['index']==ca[2]-4)]['z'].values
            ax.plot(x1,y1,z1, "o", color=colors[i%50], ms=4.0,markeredgewidth=1, markeredgecolor='black')
            ax.plot(x2,y2,z2, "o", color=colors[i%50], ms=4.0,markeredgewidth=1, markeredgecolor='black')

ls = ['' for i in range(Nc)]
fname = 'cluster_index_' +index+ '.dat'
for i in range(Nc):
    for ca in c_atms:
        if ca[0] == i+1:
            ls[i] += '{0:6d} {1:6d}'.format(ca[1]-2, ca[2]-3)
            try:
                chain.remove(ca[1]-2)
                chain.remove(ca[2]-3)
                #print(ca[1], ca[2],chain)
            except:
                pass
    with open(fname, 'at') as f:
        f.write('[ cluster No.{0:3d} ] '.format(i+1) + ls[i]+'\n')
print('unclustered')
print(chain)

fname = 'unclustered_' + ifname.split('_')[1]
il = 0
for c in chain:
    il += 1 
with open(fname, 'wt') as f:
    f.write('unclustered polymer chains\n')
    f.write('{0:8d}\n'.format(il))
for c in chain:
    nmol = int(c/1405) + 1
    c0 = c%(100000-1)
    print(c,c0)
    l = '{0:5d}PVA20   c3{1:5d}'.format(nmol, c0)
    print(l)
    x = backbones[(backbones['index']==c-1)]['x'].values[0]
    y = backbones[(backbones['index']==c-1)]['y'].values[0]
    z = backbones[(backbones['index']==c-1)]['z'].values[0]
    lr = '{0:8.3f}{1:8.3f}{2:8.3f}\n'.format(x,y,z)
    with open(fname, 'a+') as f:
        f.write(l+lr)
with open(fname, 'a+') as f:
    f.write('{0:7.5f}\t{0:7.5f}\t{0:7.5f}'.format(box))

for ca in c_atms:
    backbones.drop(backbones.index[backbones.serial==ca[1]-2], inplace=True)
    backbones.drop(backbones.index[backbones.serial==ca[2]-3], inplace=True)


xs = ['' for i in range(Nchain)]
ys = ['' for i in range(Nchain)]
zs = ['' for i in range(Nchain)]
for i in range(1, Nchain+1):
    #print(backbones[(backbones['resSeq']==i)])
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
