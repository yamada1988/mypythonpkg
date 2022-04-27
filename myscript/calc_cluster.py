import numpy as np
import sys
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import mdtraj as md
from sklearn.cluster import KMeans

colordict = matplotlib.colors.cnames
colors =  list(colordict.keys())
print(colors)


fname = 'hbonds_chain_long_others.dat'
box = 16.78288
r0 = 0.70e0
nbonds_min = 10
Nmax = 8000

with open(fname, 'rt') as f:
    hbonds = [[int(line.split()[0]), int(line.split()[1]), int(line.split()[2]), int(line.split()[3]), int(line.split()[4]), 
               float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]

N_hbond = len(hbonds)

l_ = ''
il = 0
for h in hbonds:
    l = '{0:5d} {1:5d} '.format(h[3], h[4])
    l_ += l
    il += 1
    if il == 12:
        l_ += '\n'
        with open('hbonds.ndx', 'a+') as f:
            f.write(l_)
        l_ = ''
        il = 0   

import networkx as nx
G = nx.Graph()
Nchain = 50
for i in range(1,Nchain+1):
    G.add_node(i)

for h in hbonds:                
    G.add_edge(h[1], h[2])

pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
pyplot.show() 

#sys.exit()

clusters = [[[],[]]]
id_clusts = []
id_clust = 0
lists = list(range(N_hbond+1))
for id_ in range(Nmax+1):
    print('index:', id_+1)
    index = id_+1
    try:
        chain_index1 = hbonds[id_][1]
        chain_index2 = hbonds[id_][2]
    except:
        break
    r1 = np.array(hbonds[id_][5:8])
    for l in lists:
        r2 = np.array(hbonds[l-1][5:8])
        d = np.linalg.norm(r1-r2)
        if 0.0 < d < r0:
            chain_index3 = hbonds[l-1][1]
            chain_index4 = hbonds[l-1][2]
            print(l, chain_index1, chain_index2, chain_index3, chain_index4, d)
            clusters[id_clust][0].append(id_+1)
            clusters[id_clust][1].append(chain_index1)
            clusters[id_clust][1].append(chain_index2)
            clusters[id_clust][0].append(l)
            clusters[id_clust][1].append(chain_index3)
            clusters[id_clust][1].append(chain_index4)
           
    # remove deprecate index in clusters and cluster_chain
    clusters[id_clust][0] = list(set(clusters[id_clust][0]))
    clusters[id_clust][1] = list(set(clusters[id_clust][1]))

    if clusters[id_clust][0] == []:
        continue
    if id_clust == len(id_clusts):
        clusters.append([[],[]])
    id_clusts.append(id_clust)
    id_clusts = list(set(id_clusts))
    id_clust = len(id_clusts)

clusters.remove([[],[]])
print(clusters)
print(id_clusts)

outf = 'old_hbonds_01.dat'
with open(outf, 'wt') as f:
    f.write('#clust\thbond\tchain1\tchain2\tx\ty\tz\n')
id_c = 1
for clust in clusters:
    for id_h in clust[0]:
        with open(outf, 'a+') as f:
            l = '{0:3d}\t{1:5d}\t{2:4d}\t{3:4d}\t{4:7.4f}\t{5:7.4f}\t{6:7.4f}\n'.format(id_c, id_h, hbonds[id_h-1][1], hbonds[id_h-1][2],hbonds[id_h-1][5], hbonds[id_h-1][6], hbonds[id_h-1][7])
            f.write(l)
    id_c += 1


new_clusters = [['',''] for i in range(id_clusts[-1]+1)]
dels= []
for id_c1 in range(id_clusts[-1]+1):
    print('=====cluster index:', id_c1+1)
    k = 0
    while True:
        c1 = clusters[id_c1][0]
        N0 = len(c1)
        print('=====iteration index:', k+1)
        if id_c1 in dels:
            print('This cluster already appended.')
            new_clusters[id_c1][0] = ''
            break
        c1_chain = clusters[id_c1][1]
        for id_c2 in range(id_clusts[-1]+1):
            c2 = clusters[id_c2][0]
            c2_chain = clusters[id_c2][1]
            if c1 == c2:
                continue
            #print(c1,c2)
            if list(set(c2) & set(c1)) != []:
                c1.extend(c2)
                c1_chain.extend(c2_chain)
                dels.append(id_c2)
        new_clusters[id_c1][0] = list(set(c1))
        clusters[id_c1][0] = new_clusters[id_c1][0]
        new_clusters[id_c1][1] = list(set(c1_chain))
        clusters[id_c1][1] = new_clusters[id_c1][1]
        print(N0, len(clusters[id_c1][0]))
        if len(clusters[id_c1][0]) == N0:
            break
        k+=1

while True:
    try:
        new_clusters.remove(['',''])
    except:
        break


print('before')
for nc in new_clusters:
    print(nc)
def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]
new_clusters = get_unique_list(new_clusters)
print('after')
for nv in new_clusters:
    print(nc)


outf = 'old_hbonds_02.dat'
with open(outf, 'wt') as f:
    f.write('#clust\thbond\tchain1\tchain2\tx\ty\tz\n')
for id_c, clust in enumerate(new_clusters):
    for id_h in clust[0]:
        with open(outf, 'a+') as f:
            l = '{0:3d}\t{1:5d}\t{2:4d}\t{3:4d}\t{4:7.4f}\t{5:7.4f}\t{6:7.4f}\n'.format(id_c+1, id_h, hbonds[id_h-1][1], hbonds[id_h-1][2],hbonds[id_h-1][5], hbonds[id_h-1][6], hbonds[id_h-1][7])
            f.write(l)


for j in range(10):
    for i in range(len(new_clusters)):
        try:
            if new_clusters[i][0] == '' or len(new_clusters[i][0]) < nbonds_min:
                del new_clusters[i]
        except:
            break

for i in range(len(new_clusters)):
    print(i+1, new_clusters[i][1], new_clusters[i][0])

outf = 'cluster_hbonds.dat'
with open(outf, 'wt') as f:
    f.write('#clust\thbond\tchain1\tchain2\tatm1\tatm2\tx\ty\tz\n')
id_c = 1
for clust in new_clusters:
    for id_h in clust[0]:
        with open(outf, 'a+') as f:
            l = '{0:3d}\t{1:5d}\t{2:4d}\t{3:4d}\t{4:5d}\t{5:5d}\t{6:7.4f}\t{7:7.4f}\t{8:7.4f}\n'.format(id_c, id_h, hbonds[id_h-1][1], hbonds[id_h-1][2], hbonds[id_h-1][3], hbonds[id_h-1][4],
                                                                                        hbonds[id_h-1][5], hbonds[id_h-1][6], hbonds[id_h-1][7])
            f.write(l)
    id_c += 1


Nc = len(new_clusters)
linked_nodes = [0 for i in range(Nc)]
for ic1, c1 in enumerate(new_clusters):
    for ic2 in range(ic1+1, Nc):
        if ic2 >= Nc-1:
            break
        print(ic1+1,ic2+1)
        c2 = new_clusters[ic2+1]
        if list(set(c2[1]) & set(c1[1])) != []:
            print(ic1+1, ic2+1, 'bonded')
            linked_nodes[ic1] += 1
            linked_nodes[ic2] += 1

outf = 'clusters.dat'
with open(outf, 'wt') as f:
    f.write('#clust\tbonds\txsi\tk\tlinks\tx\ty\tz\n')

pos_clust = [[0.0e0, 0.0e0, 0.0e0] for i in range(len(new_clusters))]
hbond_ndx = ['' for i in range(len(new_clusters))]
with open(outf, 'at') as f:
    for ic, clust in enumerate(new_clusters):
        N = len(clust[0])
        il = 0
        for id_h in clust[0]:
            pos_clust[ic][0] += hbonds[id_h-1][5]
            pos_clust[ic][1] += hbonds[id_h-1][6]
            pos_clust[ic][2] += hbonds[id_h-1][7]
            l = '{0:5d} {1:5d} '.format(hbonds[id_h-1][3], hbonds[id_h-1][4])
            hbond_ndx[ic] += l
            il += 1
            if il == 12:
                hbond_ndx[ic] += '\n'
                il = 0
        pos_clust[ic][0] /= N
        pos_clust[ic][1] /= N
        pos_clust[ic][2] /= N
        r = 0.0e0
        for id_h in clust[0]:
            for i in range(5,8):
                r0 = abs(hbonds[id_h-1][i] - pos_clust[ic][i-5])
                r = max([r, r0]) 
        r *= 2.0e0
        f.write('{0:4d}\t{1:3d}\t{2:5.3f}\t{3:2d}\t{4:3d}\t{5:7.4f}\t{6:7.4f}\t{7:7.4f}\n'.format(ic+1, len(clust[0]), r, len(clust[1]), linked_nodes[ic], 
                pos_clust[ic][0], pos_clust[ic][1], pos_clust[ic][2]))


with open('hbonds_clust.ndx', 'wt') as f:
    for ic in range(len(new_clusters)):
        f.write('[ cluster {0:3d} ]\n'.format(ic))
        f.write(hbond_ndx[ic])
        f.write('\n')



# Plot backbones
fig = pyplot.figure(figsize=(12,8))
ax = fig.gca(projection='3d')

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-0.0, 1.0*box)
ax.set_ylim(-0.0, 1.0*box)
ax.set_zlim(-0.0, 1.0*box)


for i in range(Nc):
    x = pos_clust[i][0]
    y = pos_clust[i][1]
    z = pos_clust[i][2]
    print('index:',i, x, y, z)
    ax.plot([x], [y], [z], "o", color=colors[i%50], ms=25)
    ax.text(x-0.20, y-0.20, z-0.20, '{0:2d}'.format(i+1), color='white' )

for ic1, c1 in enumerate(new_clusters):
    for ic2 in range(ic1+1, Nc+1):
        x = []
        y = []
        z = []
        if ic2 >= Nc-1:
            break 
        print(ic1+1, ic2+1)
        c2 = new_clusters[ic2+1]
        if list(set(c2[1]) & set(c1[1])) != []:
            print(ic1+1, ic2+1, 'bonded')
            print(pos_clust[ic1], pos_clust[ic2])
            a = np.array(pos_clust[ic1])
            b = np.array(pos_clust[ic2])
            d = np.linalg.norm(a-b)
            if d > box * 0.50e0:
                continue
            x.append(pos_clust[ic1][0])
            x.append(pos_clust[ic2][0])
            y.append(pos_clust[ic1][1])
            y.append(pos_clust[ic2][1])
            z.append(pos_clust[ic1][2])
            z.append(pos_clust[ic2][2]) 
            ax.plot(x,y,z, color='k', linestyle=':', linewidth=1)
            
pyplot.savefig('clust.eps')
pyplot.show()
