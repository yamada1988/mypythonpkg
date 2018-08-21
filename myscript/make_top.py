import shutil
import sys

args = sys.argv
inpf = args[1]
num = int(args[2])

with open(inpf, 'rt') as f:
    lines = [line.strip() for line in f]

zeroindex = lines.index('[ atoms ]')
finalindex = lines.index('[ bonds ]')

print(zeroindex)
print(finalindex)

count = 0
for i in range(zeroindex+2, finalindex):
    if lines[i]:
        count +=1
        if count >= num:
            l_ = lines[i].split()
            l_[1] += '_'
            lines[i] = '\t'.join(l_)

shutil.copyfile(inpf, inpf+'.bak')

testf = inpf
with open(testf, 'wt') as f:
    for l in lines:
        f.write(l+'\n')
