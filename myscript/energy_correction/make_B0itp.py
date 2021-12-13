import math
import shutil

#
# This script check total charge in your .itp file
# It also create "new.itp" file whose total charge equals zero.
#

inpf = 'PE100.itp'

f = open(inpf, 'rt')
lines = [line.strip() for line in f]
f.close()

zeroindex = lines.index('[ moleculetype ]')
tindex  = zeroindex + 2
lines[tindex] = 'B0' + lines[tindex]

zeroindex = lines.index('[ atoms ]')
finalindex = lines.index('[ bonds ]')

print(zeroindex)
print(finalindex)

charges = []
names = []
molnames = []
for i in range(zeroindex+2, finalindex):
    if lines[i]:
        charges.append(float(lines[i].split()[6]))
        names.append(lines[i].split()[1])
        molnames.append(lines[i].split()[3])


targetindex = 131

for i in range(targetindex, finalindex):
    if lines[i]:
        j = i - (zeroindex+2)
        index = lines[i].split()[0]
        moltype = lines[i].split()[2]
        line_sp1 = lines[i].split()[4:6]
        print(index, line_sp1)
        name = names[j] + '0'
        molname = 'B0' + molnames[j]
        charge = 0.0*charges[j]
        line_sp2 = [index, name, moltype, molname]
        line_sp2.extend(line_sp1)
        line_sp2.append(str(charge))
        print(line_sp2)
        lines[i] = '\t'.join(line_sp2)
        print(lines[i])

shutil.copyfile(inpf, inpf+'.bak')

testf = 'new.itp'
with open(testf, 'wt') as f:
    for i in range(len(lines)):
        f.write(lines[i]+'\n')
