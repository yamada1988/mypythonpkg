import numpy as np
import sys
import math

eneses = []
xs = []
aves = []
stdevs = []
errs = []

for i in range(20, 70+1, 5):
    outf = 'ave_{0:03d}.dat'.format(i)
    with open(outf, 'rt') as f:
        enes = [float(line.split()[1]) for line in f]
    x = [i]*len(enes)
    xs.append(x)
    eneses.append(enes)
    es = np.array(enes)
    #print(es)

    ave = np.mean(es)
    stdev = np.std(es)
    err = 2.0e0*stdev/math.sqrt(len(es))
    aves.append(ave)
    stdevs.append(stdev)
    errs.append(err)
    print('mean:', ave)
    print('stdev:', stdev)
    print('95%err:', err)

outf = 'aveuv.tt'
with open(outf, 'wt') as f:
    for ix, x in enumerate(xs):
        print(ix,x)
        for ex in eneses[ix]:
            f.write('{0:d}\t{1:8.5f}\n'.format(x[0], ex))
        
outf2 = 'ave_aveuv.tt'
with open(outf2, 'wt') as f:
    for ix, x in enumerate(xs):
        f.write('{0:d}\t{1:8.5f}\t{2:6.4f}\n'.format(x[0], aves[ix],errs[ix]))

outf3 = 'averaged_ave.dat'
with open(outf3, 'rt') as f:
    eneses = [float(line.split()[1]) for line in f]
eneses = np.array(eneses)
ave_aveuv = np.mean(eneses)
print(ave_aveuv)

outf4 = 'TotalAverage.dat'
with open(outf4, 'wt') as f:
    f.write('{0:7.5f}\n'.format(ave_aveuv))

print(np.std(aves))
