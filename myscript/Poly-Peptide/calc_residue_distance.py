import mdtraj as md
import numpy as np
import itertools
import sys

args = sys.argv
index = int(args[1])

fname = "MD/chain{0:04d}.xtc".format(index)
#fname = "MD/test{0:04d}.xtc".format(index)
sysgro = "SYS/chain{0:04d}.gro".format(index)
ts = md.iterload(fname,top=sysgro)
for i,t_test in enumerate(ts):
    if i == 0:
        t = t_test
    print(i)

N_iterloads = i+1

top = t.topology
N = 48
atms = []
for i in range(N):
    line = "resSeq {0:d} and name C".format(i+1)
    atms.append(top.select(line)[0])

pairs = np.array((list(itertools.permutations(atms, 2))))
#print(pairs)

distance_matrix = np.zeros((N,N-1))
ts = md.iterload(fname,top=sysgro) #reinitializeation
for i,t in enumerate(ts):
    #print(i,t)
    print(i)
    dm = md.compute_distances(t, pairs)
    #print(dm)
    dm = np.mean(dm, axis=0)
    #print(dm)
    dm = dm.reshape(N,N-1)
    distance_matrix += dm
    #if i == 2:
    #    break

distance_matrix /= N_iterloads
#print(distance_matrix)

# insert diagonal element (0.0e0) in distance_matrix
output_distance_matrix = []

for i,d in enumerate(distance_matrix):
    #print(i, d)
    output_distance_matrix.append(np.insert(d, i, 0.0e0))

output_distance_matrix = np.array(output_distance_matrix)
#print(output_distance_matrix)  


ofname = "DAT/residue_distance{0:04d}.dat".format(index)
np.savetxt(ofname, output_distance_matrix, fmt="%.4e")
