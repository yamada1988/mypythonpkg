import math
import shutil

#
# This script check total charge in your .itp file
# It also create "new.itp" file whose total charge equals zero.
#

inpf = 'PVP100.itp'

f = open(inpf, 'rt')
lines = [line.strip() for line in f]
f.close()

zeroindex = lines.index('[ atoms ]')
finalindex = lines.index('[ bonds ]')

print(zeroindex)
print(finalindex)

charges = []
names = []
for i in range(zeroindex+2, finalindex):
    if lines[i]:
        charges.append(float(lines[i].split()[6]))
        names.append(lines[i].split()[1])


for i in range(zeroindex+2, finalindex):
    if lines[i]:
        j = i - (zeroindex+2)
        index = lines[i].split()[0]
        line_sp1 = lines[i].split()[2:6]
        print(index, line_sp1)
        name = names[j] + '0'
        charge = 0.0*charges[j]
        line_sp2 = [index, name]
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
