import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans

fname = 'hbonds_chain.dat'
num = 48
box = 16.50
colordict = matplotlib.colors.cnames
colors =  list(colordict.keys())
print(colors)


with open(fname, 'rt') as f:
    crystalines = [[float(line.split()[5]), float(line.split()[6]), float(line.split()[7])] for line in f if not line.startswith('#')]
    f.seek(0)
    crystalines_02 = [[float(line.split()[8]), float(line.split()[9]), float(line.split()[10])] for line in f if not line.startswith('#')]

crystalines.extend(crystalines_02)
crystalines = np.array(crystalines)
pred = KMeans(n_clusters=num).fit_predict(crystalines)
print(pred)

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
    ax.scatter(c[0], c[1], c[2], marker="*", color=colors[pred[i]], s=10)

plt.show()
