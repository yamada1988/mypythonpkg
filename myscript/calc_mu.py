import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import os
import sys
import math
sns.set(font='Kozuka Gothic Pro', style="whitegrid")
#sns.color_palette("hls",24)
sns.set_palette("husl", 12)

#
# Run sh/make_slvfes.sh befor run this script
#

def print_ex(strname, num):
    print("{0}: {1:8.5f}".format(strname, num))

def describe_list(list):
    print_ex("    mean", np.mean(list))
    print_ex("    var", np.var(list))
    print_ex("    stdev", np.std(list))
#    print_ex("    skew", scipy.stats.skew(list))
#    print_ex("    kurt", scipy.stats.kurtosis(list))
    print_ex("    max", np.max(list))
    print_ex("    med", np.median(list))
    print_ex("    min", np.min(list))

def calc_cumave(arrayname, linenum, rawnum, fname):
    maxnum = rawnum * (linenum-1) + rawnum
    datalist = [0]*maxnum
    
    for i in range(linenum):        
        for j in range(rawnum):
            indx = rawnum * i + j
            datalist[indx] = arrayname[i][j]


    with open(fname, 'wt') as f:
        f.write('# mu (kcal/mol)\terror\n')

    for indx in range(maxnum):
        ave_indx = np.mean(datalist[:indx+1])
        std_indx = np.std(datalist[:indx+1])
        with open(fname, 'a+') as f:
            f.write('{0:d}\t{1:7.5f}\t{2:7.5f}\n'.format(indx, ave_indx, std_indx*2/math.sqrt(indx)))
hists = []
Nmax = 50
Nmin = 39
indexes = range(Nmax, Nmin-1, -1)
Ntot = len(indexes)
for indx in indexes:
    filename = "MELT_{0:03d}/ERmods/log.tt".format(indx)
    figname = "DAT/mu{0:03d}.eps".format(indx)
    if not os.path.exists(filename):
        continue
    else:
        with open(filename, 'rt') as f:
            mus = [map(float, line[5:].split()) for line in f]

    mus = np.array(mus)
    rawnum = len(mus[0])
    linenum = len(mus)
    calc_cumave(mus, linenum, rawnum, 'MELT_{0:03d}/ERmods/slvfes.tt'.format(indx))

    print("index:{0:03d}".format(indx))
    print("""========================= 
  Statistics of Energies 
=========================""")
    describe_list(mus)

    hist, bins = np.histogram(mus, range=(-8.0e0, 1.0e0), bins=100, normed=True)
    hists.append(hist)

for i in range(Ntot):
    plt.step(bins[:-1], hists[i], where='post', label='{0:03d}'.format(Nmax-i))
#    plt.plot(bins[:-1], hists[i], label='{0:03d}'.format(Nmax-i))
    with open('ermod_010conf_0500ps_rand_40000frame/slvfe_{0:03d}.dat'.format(Nmax-i), 'w+') as f:
        for j in range(len(bins)-1):
            f.write('{0:11.5f}\t{1:18.17f}\n'.format(bins[j], hists[i][j]))

plt.legend(ncol=3)
plt.show()
  
