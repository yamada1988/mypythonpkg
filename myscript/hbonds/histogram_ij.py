import numpy as np
import pickle
import sys

def load_dumps(f):
    obj = []
    while 1:
        try:
            obj.append(pickle.load(f))
        except:
            break
    return obj

inpdir = 'sij_md/'
with open(inpdir + 'sij.pickle', 'rb') as f:
    data = load_dumps(f)

del data[0]
for i,d  in enumerate(data):
    print(i, d)


T = len(data)
t=0
hist = np.array([0.0 for d in data])
for i in range(5000):
    jlist = list(set([x for d in data for x in list(set(d[1][i].indices))]))
    print("###\n", i, jlist)
    for j in jlist:
        print(j)
        ijpairs = [d[1][i,j] for d in data]
        for ijp in ijpairs:
            if ijp == 1.0:
                t+=1
            elif ijp != 1:
                if t != 0:
                    hist[t] += 1
                    print(t)
                    t=0

hist /= np.sum(hist)
with open('hist.dat', 'wt') as f:
    for it in range(1,T):
        #print(data[it][0], p[it])
        l = '{0:7.3f}\t{1:6.5e}\n'.format(data[it][0], hist[it])
        f.write(l)
