import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.animation as animation
import seaborn as sns
import mdtraj as md
import numpy as np
import sys

# Parameters

Nmol = 1611
Nchain = 50
Natm = 705
Nseg = 10
Nstart = 4654

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

sns.set_palette("hls", Nchain)
hbondf = 'hbonds_par{0:04d}_chain{1:05d}.dat'.format(index_, frame_)

# Hydrogen bonding between inter-chain
with open(hbondf, 'rt') as f:
    atm_bonded = [int(line.split()[3])-3 for line in f if not line.startswith('#')]


for a in atm_bonded:
    print(a)
    chn_index = int(a/Natm)+1
    atm_index = a%Natm
    print(chn_index, atm_index)
    if atm_index <= 71:
        a_index = 1
    elif atm_index >= 632:
        a_index = Nseg
    else:
        a_index = int(atm_index/71)
    g_index = (chn_index-1)*Nseg+a_index
    print(a_index, g_index) 

G = nx.Graph()
for i in range(1,Nchain*Nseg+1):
    G.add_node(i)

for j in range(1, Nchain+1):
    pair_bonded = [(i,i+1) for i in range((j-1)*Nseg+1,j*Nseg)]
    print(pair_bonded)
    G.add_edges_from(pair_bonded)

nx.draw_networkx(G)
plt.show()
