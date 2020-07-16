import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

from mpl_toolkits.mplot3d import Axes3D

Nt = 1500
for it in range(Nt, Nt+1,50):
    fname = 'DAT/angle_ave/rho_hr{0:05d}.dat'.format(it)
    print(fname)
    rs = np.loadtxt(fname)
    rs = rs[rs[:,3] > 0.960]

    X = rs[::4,0]
    Y = rs[::4,1]
    Z = rs[::4,2]
    value = rs[::4,3]

    fig = plt.figure()
    ax = Axes3D(fig)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    cm = plt.cm.get_cmap('bwr')
    mappable = ax.scatter(X, Y, Z, c=value, cmap=cm, s=40, vmin=-0.6,vmax=1.3)
    fig.colorbar(mappable, ax=ax)

plt.show()
