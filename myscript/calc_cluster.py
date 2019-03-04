import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans

fname = 'hbonds_chain.dat'
num = 32
box = 18.915
colordict = matplotlib.colors.cnames
colors =  list(colordict.keys())
print(colors)


with open(fname, 'rt') as f:
    crystalines = [[float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]
    f.seek(0)
    crystalines_02 = [[float(line.split()[8]), float(line.split()[9]), float(line.split()[10])] for line in f if not line.startswith('#')]

crystalines.extend(crystalines_02)
crystalines = np.array(crystalines)
cls = KMeans(n_clusters=num)
pred = cls.fit_predict(crystalines)
print(pred)
centers = cls.cluster_centers_

for i in range(num):
    print(i, list(pred).count(i))

 # Plot backbones
fig = pyplot.figure(figsize=(12,8))
ax = Axes3D(fig)

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

ax.set_xlim(-0.0, 1.0*box)
ax.set_ylim(-0.0, 1.0*box)
ax.set_zlim(-0.0, 1.0*box)

# Plot hydrogen bonding
#print(x_ac_bonded)
for i,c in enumerate(crystalines):
    #print(i,c)
    ax.scatter(c[0], c[1], c[2], marker='*', c=colors[pred[i]], edgecolors=colors[pred[i]], s=50)

outf = 'center_bonds.dat'
with open(outf, 'wt') as f:
    f.write('# index\tx\ty\tz\n')
    for j in range(num):
        ax.scatter(centers[j, 0], centers[j, 1], centers[j,2], marker='*', c=colors[j], s=500, edgecolors='black')
        f.write('{0:3d}\t{1:5.3f}\t{2:5.3f}\t{3:5.3f}\n'.format(j, centers[j,0], centers[j,1], centers[j,2]))

plt.show()
