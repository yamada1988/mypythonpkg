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
#for i,d  in enumerate(data):
    #print(i, d)

def calc_hist(i, t):
    #print(i)
    jlist = list(set([x for x in list(set(data[t][1][i].indices))]))
    jlist.sort()
    if jlist != []:
        histo_numpair[i] += jlist
        histo_numpair[i] = list(set(histo_numpair[i]))
    #print(jlist)
    #print(hist_numpair[len(jlist)], numid[i])
    #numid[i] = len(jlist)
    return jlist


T = len(data)
import os
try:
    os.mkdir('travels_num')
except:
    pass


histo_numpair = [[] for i in range(1000)]
alpha = 0.10

nmin = anim_index
for i in range(nmin, nmin+1):
    with open('travels_num/water{0:04d}.dat'.format(i+1), 'wt') as f:
        f.write('#time(ns)\tnum\tpairs\n')
    for t in range(T):
        numid = calc_hist(i, t)
        if numid == []:
            ls = '{0:11.4f}'.format(t*alpha/1000.0)+'\n'
        else:
            ls = '{0:11.4f}\t'.format(t*alpha/1000.0)+'{0:>3d}'.format(len(histo_numpair[i]))
            lslen = len(ls)
            addlen = 45 - lslen
            numid = np.array(numid)
            numid += 1 # modify the index of molecules
            ls += '\t' + ' '.join(map(str, numid)) + '\n'

        with open('travels_num/water{0:04d}.dat'.format(i+1), 'a+') as f:
            f.write(ls)
