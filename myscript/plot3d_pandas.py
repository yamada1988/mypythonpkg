from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import mdtraj as md
import numpy as np
import sys

sns.set_palette("hls", 50)

index_ = int(sys.argv[1])

fname = 'npt_par{0:04d}_nopbc.gro'.format(index_)
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

print(df)

tmp = df[df['resName']=='PVA0c']
backbones = tmp[tmp['name']=='c3']

print(backbones)
fig = pyplot.figure()
ax = Axes3D(fig)

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-3.0, 1.40*box)
ax.set_ylim(-3.0, 1.40*box)
ax.set_zlim(-3.0, 1.40*box)

xs = ['' for i in range(50)]
ys = ['' for i in range(50)]
zs = ['' for i in range(50)]

for i in range(1, 50+1):
    xs[i-1] = backbones[(backbones['resSeq']==i)]['x'].values
    ys[i-1] = backbones[(backbones['resSeq']==i)]['y'].values
    zs[i-1] = backbones[(backbones['resSeq']==i)]['z'].values
    ax.plot(xs[i-1], ys[i-1], zs[i-1], "o-", ms=1)


pyplot.show()
