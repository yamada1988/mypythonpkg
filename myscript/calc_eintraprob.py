import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import os
import sys
sns.set(font='Kozuka Gothic Pro', style="whitegrid")
#sns.color_palette("hls",24)
sns.set_palette("husl", 50)

def print_ex(strname, num):
    print("{0}: {1:8.5f}".format(strname, num))

def describe_list(list):
    print_ex("    mean", np.mean(list))
    print_ex("    var", np.var(list))
    print_ex("    stdev", np.std(list))
    print_ex("    skew", scipy.stats.skew(list))
    print_ex("    kurt", scipy.stats.kurtosis(list))
    print_ex("    max", np.max(list))
    print_ex("    med", np.median(list))
    print_ex("    min", np.min(list))


e_intras = []
hists = []
cumls = []
numbins = 500
estart = 1.050e+04
eend = 1.250e+04
dx = (eend - estart) / float(numbins)

Nmax = 50
Nmin = 15
indexes = range(Nmax, Nmin-1, -1)
Ntot = len(indexes)
print(Ntot)
for indx in indexes:
    filename = "solute/eintra_{0:03d}.log".format(indx)
    figname = "DAT/e_intra{0:03d}.eps".format(indx)
    if not os.path.exists(filename):
        continue
    else:
        with open(filename, 'rt') as f:
            e_intras = [float(line.split()[4]) for line in f]
#    print(e_intras)

    e_intras = np.array(e_intras)
    print("index:{0:03d}".format(indx))
    print("""========================= 
  Statistics of Energies 
=========================""")
    describe_list(e_intras)

    hist, bins = np.histogram(e_intras, range=(estart, eend), bins=numbins, normed=True)
    cuml = np.cumsum(hist)*dx
    hists.append(hist)
    cumls.append(cuml)

for i in range(Ntot):
    plt.plot(bins[:-1], hists[i], label='{0:03d}'.format(Nmax-i))
    with open('eintra_{0:03d}.dat'.format(Nmax-i), 'w') as f:
        for j in range(len(bins)-1):
            f.write('{0:11.5f}\t{1:18.17f}\t{2:18.17f}\n'.format(bins[j], hists[i][j], cumls[i][j]))

plt.legend(ncol=3)
plt.title(r'P($\phi)')
plt.show()
plt.savefig(figname)
