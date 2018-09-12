import shutil
import sys
import os

args = sys.argv

Nmin = int(args[1])
Nmax = int(args[2])

Nblock = 5

Nst = int(args[3]) 
Ned = int(args[4])
Ni = int(args[5])

Nconf = Ned - Nst + 1
inpf_0 = 'MELT_{0:03d}/ERmod_0001/soln/engsln.01'.format(Nmin)
with open(inpf_0, 'rt') as f:
    ene_coord = [line.split()[0] for line in f]

num_lines = len(ene_coord)
for j in range(Nmax, Nmin-1, -1):
    inpdir = 'MELT_{0:03d}/'.format(j)
    esln = [[0 for i in range(Nblock+1)] for j in range(Nconf)]
    esln_ave = [ 0 for i in range(Nblock+1)]
    for il,l in enumerate(range(Nst, Ned+1)):
        for k in range(1, Nblock+1):
            inpf_sln = inpdir + 'ERmod_{0:04d}/soln/engsln.{1:02d}'.format(l,k)
            print(inpf_sln)
            with open(inpf_sln) as f:
                esln[il][k-1] = [line.split() for line in f]
#                print(esln[l-1][k-1])
                esln_ave[k-1] = ['' for line in range(len(ene_coord))]

#    print(esln[0])
    sum_ = [[0 for line in range(num_lines)] for line in range(Nblock)]
    ave = [[0 for line in range(num_lines)] for line in range(Nblock)]
    print('==========\n\n==========')
    for ie in range(num_lines):
        for k in range(1, Nblock+1):
            sum_[k-1][ie] = 0.0e0
            for il,l in enumerate(range(Nst, Ned+1)):
#                print('l:', l, 'k:', k)
#                print(esln[il][k-1][ie][2])
                sum_[k-1][ie] += float(esln[il][k-1][ie][2]) 
                ave[k-1][ie] = sum_[k-1][ie]/ Nconf
                
            esln_ave[k-1][ie] = [ene_coord[ie], '1', ave[k-1][ie]]
#            print(esln_ave[k-1][ie])
    for k in range(1, Nblock+1):
        outd0 = inpdir + 'aveERmod_{0:02d}/'.format(Ni)
        outd1 = outd0 + 'soln/'
        outf = outd1 + 'engsln.{0:02d}'.format(k)
        if not os.path.isdir(outd0):
            os.mkdir(outd0)     
        if not os.path.isdir(outd1):
            os.mkdir(outd1) 
        if os.path.isfile(outf):
            shutil.move(outf, outf+'.bak')
        with open(outf, 'a+') as f:
            for ie in range(num_lines):
                f.write(esln_ave[k-1][ie][0]+'\t'+esln_ave[k-1][ie][1]+'\t'+'{0:12.9E}'.format(esln_ave[k-1][ie][2])+'\n')

    w_slns =[0.0e0 for line in range(Nconf)]
    w_sln = [0.0e0 for line in range(Nblock)]
    for il,l in enumerate(range(Nst, Ned+1)):
        w_slnf = inpdir + 'ERmod_{0:04d}/soln/weight_soln'.format(l)
        with open(w_slnf, 'rt') as f:
            w_slns[il] = [float(line.split()[1]) for line in f]
        for k in range(1, Nblock+1):
            w_sln[k-1] += w_slns[il][k-1]
    w_sln = list(map(lambda x: x/Nconf, w_sln))
   # print(w_sln)
    outf = outd1 + 'weight_soln'
    with open(outf, 'wt') as f:
        for k in range(1, Nblock+1):
            str_ = '{0:d}\t{1:18.9E}\n'.format(k, w_sln[k-1])
            f.write(str_)

    aveuvs = [0.0e0 for line in range(Nconf)]
    aveuv = [0.0e0 for line in range(Nblock)]
    for il,l in enumerate(range(Nst, Ned+1)):
        aveuvf = inpdir + 'ERmod_{0:04d}/soln/aveuv.tt'.format(l)
        with open(aveuvf, 'rt') as f:
            aveuvs[il] = [float(line.split()[1]) for line in f]
        for k in range(1, Nblock+1):
            aveuv[k-1] += aveuvs[il][k-1]
    aveuv = list(map(lambda x: x/Nconf, aveuv))
    
    outf = outd1 + 'aveuv.tt'
    with open(outf, 'wt') as f:
        for k in range(1, Nblock+1):
            str_ = '{0:d}\t{1:8.5f}\n'.format(k, aveuv[k-1])
            f.write(str_)
