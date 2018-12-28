from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import seaborn as sns
import mdtraj as md
import numpy as np
import sys

# Parameters

Nmol = 1611
Nchain = 50
Nstart = 4654
d_hbond = 0.240 # (nm)

# Input Setting
index_ = int(sys.argv[1]) 
frame_ = int(sys.argv[2])
print('index:',index_)

fname = 'npt_par{0:04d}_nopbc{1:05d}.gro'.format(index_, frame_)
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

#print(df)
tmp = df[df['resName']=='PVA0c']
backbones = tmp[tmp['name']=='c3']
acceptor = df[df['name'] == 'oh']
donor = df[df['name'] == 'ho']

#print(acceptor)
#print(donor)

# Output Setting
sns.set_palette("hls", Nchain)
hbondf = 'hbonds_par{0:04d}_chain{1:05d}.dat'.format(index_, frame_)

# Hydrogen bonding between inter-chain
with open(hbondf, 'rt') as f:
    pos_ac_bonded = [[float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]

x_ac_bonded = list(zip(*pos_ac_bonded)[0])
y_ac_bonded = list(zip(*pos_ac_bonded)[1]) 
z_ac_bonded = list(zip(*pos_ac_bonded)[2])
# Plot backbones
fig = pyplot.figure()
ax = Axes3D(fig)

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-3.0, 1.40*box)
ax.set_ylim(-3.0, 1.40*box)
ax.set_zlim(-3.0, 1.40*box)

xs = ['' for i in range(Nchain)]
ys = ['' for i in range(Nchain)]
zs = ['' for i in range(Nchain)]

for i in range(1, Nchain+1):
    xs[i-1] = backbones[(backbones['resSeq']==i)]['x'].values
    ys[i-1] = backbones[(backbones['resSeq']==i)]['y'].values
    zs[i-1] = backbones[(backbones['resSeq']==i)]['z'].values
    ax.plot(xs[i-1], ys[i-1], zs[i-1], "o-", ms=0.3)

# Plot hydrogen bonding
#print(x_ac_bonded)
ax.plot(x_ac_bonded, y_ac_bonded, z_ac_bonded, "*", color="k", ms=6)
ax.text2D(0.05, 0.95, "{0:4.1f} ps".format(float(frame_)), transform=ax.transAxes)

ofname = 'npt_par{0:04d}_nopbc{1:05d}.eps'.format(index_, frame_)
pyplot.savefig(ofname)

#pyplot.show()
