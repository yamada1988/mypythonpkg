import sys

txtfile = sys.argv[1]

with open(txtfile, 'rt') as f:
    ene_incrs = [float(line.split()[1]) for line in f if not line.startswith('=') and not line.startswith('\n') ]
    f.seek(0)
    errors = [float(line.split()[2]) for line in f if not line.startswith('=') and not line.startswith('\n') ]
    f.seek(0)
    index0 = [int(line.split('MELT_')[1].split('/')[0]) for line in f if line.startswith('=')]

N = len(errors)
#print(mu_incrs)
#print(index)
#print(errors)

outfile = 'ene_incr.dat'
with open(outfile, 'w') as f:
    for i,ene in enumerate(reversed(ene_incrs)) :
        f.write('{0:3d}\t{1:8.5f}\t{2:6.5f}\n'.format(i+1, ene, errors[N-i-1]))
