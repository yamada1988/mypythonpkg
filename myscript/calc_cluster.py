import matplotlib
from matplotlib import pyplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pyclustering
from pyclustering.cluster import xmeans
from pyclustering.cluster import cluster_visualizer
import numpy as np


fname = 'hbonds_chain.dat'
knum = 100
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
data = crystalines
init_center = pyclustering.cluster.xmeans.kmeans_plusplus_initializer(data, 3).initialize() # 初期値決定　今回は、初期クラスタ2です
xm = pyclustering.cluster.xmeans.xmeans(data, init_center, kmax=128, ccore=False)
xm.process() # クラスタリングします

clusters = xm.get_clusters()
centers = xm.get_centers()

outf = 'center_hbonds.dat'
ic = 1
with open(outf, 'wt') as f:
    f.write('# index\tnum\tx\ty\tz\n')
    for c in centers:
        nc = len(clusters[ic-1])
        f.write('{0:3d}\t{1:4d}\t{2:6.3f}\t{3:6.3f}\t{4:6.3f}\n'.format(ic, nc, c[0], c[1], c[2]))
        ic += 1

# Visualize clustering results
visualizer = cluster_visualizer()
visualizer.append_clusters(clusters, data)
visualizer.append_cluster(centers, None, marker='*', markersize = 100)
visualizer.show()
