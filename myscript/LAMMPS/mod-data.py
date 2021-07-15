import sys

Np = 70200
nwater = 57078
ntip3p = 3
N = Np + nwater *  ntip3p

args = sys.argv
fname = args[1]
outf = fname.split('.')[0] + '-mod.data'
with open(fname, 'rt') as f:
    lines = [line for line in f]

#print(lines)
for il,line in enumerate(lines):
    if line.startswith('Atoms'):
        atom_index_st = il + 2
    elif line.startswith('Velocities'):
        atom_index_end = il - 2
    elif line.startswith('Bonds'):
        bond_index_st = il + 2
    elif line.startswith('Angles'):
        bond_index_end = il - 2
        angle_index_st = il + 2
    elif line.startswith('Dihedrals'):
        angle_index_end = il - 2

#print(bond_index_st, bond_index_end, angle_index_st, angle_index_end)

# Atom section
atomlines = ['']*N
for il0,il in enumerate(range(atom_index_st, atom_index_end+1)):
    ls = lines[il].split()
    l0 = int(ls[0])
    if l0 <= N:
        atomlines[il0] = lines[il]
#print(atomlines)

# Bond section
Nbond = bond_index_end - bond_index_st + 1
bondlines = ['']*Nbond
ibond_mod = 0
for il0, il in enumerate(range(bond_index_st, bond_index_end+1)):
    ls = lines[il].split()
    lb = int(ls[-1])
    if lb <= N:
        bondlines[ibond_mod] = '{0:d}\t'.format(ibond_mod+1) + '\t'.join(ls[1:]) + '\n'
        ibond_mod += 1
bondlines = bondlines[:ibond_mod]
#print(bondlines)

# Angle section
Nangle = angle_index_end - angle_index_st + 1
anglelines = ['']*Nangle
iangle_mod = 0
for il0, il in enumerate(range(angle_index_st, angle_index_end+1)):
    ls = lines[il].split()
    la = int(ls[-1])
    if la <= N:
        anglelines[iangle_mod] = '{0:d}\t'.format(iangle_mod+1) + '\t'.join(ls[1:]) + '\n'
        iangle_mod += 1
anglelines = anglelines[:iangle_mod]
#print(anglelines)

# Output
lines[2] = '{0:6d} atoms\n'.format(N)
lines[4] = '{0:6d} bonds\n'.format(ibond_mod)
lines[6] = '{0:6d} angles\n'.format(iangle_mod)

with open(outf, 'wt') as f:
    for line in lines[:atom_index_st]:
        f.write(line)

    for line in atomlines:
        f.write(line)
   
    f.write('\n'+'Bonds\n\n')
    for line in bondlines:
        f.write(line)

    f.write('\n'+'Angles\n\n')
    for line in anglelines:
        f.write(line)
 
    for line in lines[angle_index_end+1:]:
        f.write(line)
