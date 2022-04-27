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
print('frame:',frame_)

fname = 'npt_par{0:04d}_nopbc{1:05d}.gro'.format(index_, frame_)
t = md.load(fname)
top = t.topology
df, b = top.to_dataframe()
pos = t.xyz
box = t.unitcell_lengths[0,0]

df['x'] = pos[0,:,0]
df['y'] = pos[0,:,1]
df['z'] = pos[0,:,2]

print(df)
tmp = df[df['resName'] == 'HOH']
sol = tmp[tmp['name'] == 'O']
x_sol_bonded = list(sol['x'])
y_sol_bonded = list(sol['y'])
z_sol_bonded = list(sol['z'])
Nhoh = len(x_sol_bonded)
print(x_sol_bonded)

tmp = df[df['resName']=='PVA1c']
backbones = tmp[tmp['name']=='c3'].reset_index()[1::2]

sns.set_palette("hls", Nchain)
hbondf = 'hbonds_par{0:04d}_chain{1:05d}.dat'.format(index_, frame_)


# Calc Ow-crosslink point distance
with open(hbondf, 'rt') as f:
    atm_bonded = [int(line.split()[3])-2 for line in f if not line.startswith('#')]

hbonded = backbones[backbones['serial'].isin(atm_bonded)]
x_ac_bonded = list(hbonded['x'])
y_ac_bonded = list(hbonded['y'])
z_ac_bonded = list(hbonded['z'])
Ncrosslink = len(x_ac_bonded)

x_gel = []
y_gel = []
z_gel = []
for i in range(Ncrosslink):
    r = np.array([x_ac_bonded[i], y_ac_bonded[i], z_ac_bonded[i]])
    print(r)
    for ih in range(Nhoh):
        r0 = np.array([x_sol_bonded[ih], y_sol_bonded[ih], z_sol_bonded[ih]])
        dr = np.abs(r - r0)
        dr -= np.round(dr/box)*box
        d = np.sqrt(np.sum(dr**2))
        if d < 0.50:
            print(i,ih)
            x_gel.append(r0[0])
            y_gel.append(r0[1])
            z_gel.append(r0[2])

# Plot backbones
fig = pyplot.figure()
ax = Axes3D(fig)

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-2.0, 1.20*box)
ax.set_ylim(-2.0, 1.20*box)
ax.set_zlim(-2.0, 1.20*box)

xs = ['' for i in range(Nchain)]
ys = ['' for i in range(Nchain)]
zs = ['' for i in range(Nchain)]

for i in range(1, Nchain+1):
    xs[i-1] = backbones[(backbones['resSeq']==i)]['x'].values
    ys[i-1] = backbones[(backbones['resSeq']==i)]['y'].values
    zs[i-1] = backbones[(backbones['resSeq']==i)]['z'].values
    ax.plot(xs[i-1], ys[i-1], zs[i-1], "o-", lw=0.20, ms=4.0)


# Plot hydrogen bonding
#print(x_ac_bonded)
ax.plot(x_ac_bonded, y_ac_bonded, z_ac_bonded, "*", color="k", ms=5)
ax.text2D(0.90, 0.95, "{0: 4.1f} ps".format(float(frame_)), transform=ax.transAxes)

ax.plot(x_gel, y_gel, z_gel, "s", color="k", ms=2)




ofname = 'npt_par{0:04d}_bead{1:05d}.eps'.format(index_, frame_)
pyplot.savefig(ofname)

pyplot.show()
