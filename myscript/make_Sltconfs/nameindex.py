import sys

args = sys.argv
fname = args[1]

with open(fname, 'rt') as f:
    lines = [line.strip() for line in f]

numatm = 477
index0 = 7
index1 = numatm + index0 - 1

resnames = []
resnums = []
resindices = []
k0 = 1
for i in range(index0, index1):
    #print(i,lines[i])
    li = lines[i].split()
    k = int(li[2])
    rname = li[3]
    print(li[3],k,k0)
    if k == k0:
        resnames.append(rname)
        resnums.append(k)
        resindices.append(int(li[0])-1)
        k0 += 1

#print(resnames, resnums, resindices)

for i0, rn in enumerate(resnames):
    print(i0-1,rn)
