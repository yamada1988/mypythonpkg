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
    print_ex("  95%err", 2.0e0*np.std(list)/math.sqrt(len(list)))
    print_ex("    var", np.var(list))
    print_ex("    stdev", np.std(list))
    print_ex("    max", np.max(list))
    print_ex("    med", np.median(list))
    print_ex("    min", np.min(list))

def calc_cumave(arrayname, linenum, fname):
    maxnum = linenum
    datalist = [0]*maxnum
    
    for i in range(linenum):        
        indx = i
        datalist[indx] = arrayname[i]


    with open(fname, 'wt') as f:
        f.write('# aveuv (kcal/mol)\terror\n')

    for indx in range(maxnum):
        ave_indx = np.mean(datalist[:indx+1])
        std_indx = np.std(datalist[:indx+1])
        with open(fname, 'a+') as f:
            f.write('{0:d}\t{1:8.5f}\t{2:8.5f}\n'.format(indx, ave_indx, std_indx*2.0e0/math.sqrt(indx)))


hists = []
Nmax = 50
Nmin = 39
indexes = range(Nmax, Nmin-1, -1)
Ntot = len(indexes)
for indx in indexes:
    mus = []
    filename = "MELT_{0:03d}/ERmods/loguv.tt".format(indx)
    figname = "DAT/aveuv{0:03d}.eps".format(indx)
    if not os.path.exists(filename):
        continue
    else:
        with open(filename, 'rt') as f:
            euvs = [float(line[9:].split()[0]) for line in f]

    euvs = np.array(euvs)
    print(euvs)
    linenum = len(euvs)

#    print(sum(euvs))  
 
    calc_cumave(euvs, linenum, 'MELT_{0:03d}/ERmods/aveuv.tt'.format(indx))

    print("index:{0:03d}".format(indx))
    print("""========================= 
  Statistics of Energies 
=========================""")
    describe_list(euvs)

    hist, bins = np.histogram(euvs, range=(-5.0e0, -2.0e0), bins=200, normed=True)
    hists.append(hist)

#for i in range(Ntot):
#    plt.step(bins[:-1], hists[i], where='post', label='{0:03d}'.format(Nmax-i))
#    plt.plot(bins[:-1], hists[i], label='{0:03d}'.format(Nmax-i))

#plt.legend(ncol=3)
#plt.show()
