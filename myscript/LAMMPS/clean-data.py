import sys

args = sys.argv
fname = args[1]
outf = fname.split('.')[0] + '-clean.data'
with open(fname, 'rt') as f:
    lines = [line for line in f]

#print(lines)
for il,line in enumerate(lines):
    if line.startswith('Atoms'):
        atom_index_st = il + 2
    elif line.startswith('Velocities'):
        atom_index_end = il - 2

print(atom_index_st, atom_index_end)
N = atom_index_end - atom_index_st + 1
ls = lines[atom_index_st].split()
l0 = ls[0]
lr = '\t'.join(ls[0:-3])+'\n'
#print(l0, lr)
d = {l0:lr}
for l in lines[atom_index_st+1:atom_index_end+1]:
    #print(l)
    ls = l.split()
    l0 = ls[0]
    lr = '\t'.join(ls[0:-3])+'\n'
    d[l0] = lr

#print(d)
i0 = atom_index_st
for i in range(N):
    #print(i)
    lines[i+i0] = d['{0:d}'.format(i+1)]
    #print(atomlines[i])


#sys.exit()

with open(outf, 'wt') as f:
    for il, line in enumerate(lines):
        f.write(line)
