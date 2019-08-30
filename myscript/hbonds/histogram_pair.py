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

inpdir = 'sij_md_skip/'
with open(inpdir + 'sij.pickle', 'rb') as f:
    data = load_dumps(f)

del data[0]
for i,d  in enumerate(data):
    print(i, d)

def calc_hist(i, t):
    print(i)
    jlist = list(set([x for x in list(set(data[t][1][i].indices))]))
    jlist.sort()
    print(jlist)
    #print(hist_numpair[len(jlist)], numid[i])
    return jlist

T = len(data)
import os
try:
    os.mkdir('travels_num')
except:
    pass


for i in range(0,1):
    with open('travels_num/water{0:04d}.dat'.format(i+1), 'wt') as f:
        f.write('#time\tpair_num\n')
    for t in range(T):
        numid = calc_hist(i, t)
        if numid == []:
            ls = '\n'
        else:
            ls = '\t'.join(map(str, numid))
            ls += '\n'

        with open('travels_num/water{0:04d}.dat'.format(i+1), 'a+') as f:
            f.write(ls)
