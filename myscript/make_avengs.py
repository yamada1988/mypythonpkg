import shutil

Nmax = 50
Nmin = 50

Nblock = 5

Nst =  1
Ned =  20
Nconf = Ned - Nst + 1
inpf_0 = 'MELT_050/ERmod_0001/soln/engsln.01'
with open(inpf_0, 'rt') as f:
    ene_coord = [line.split()[0] for line in f]

num_lines = len(ene_coord)
for j in range(Nmax, Nmin-1, -1):
    inpdir = 'MELT_{0:03d}/'.format(j)
    esln = [[0 for i in range(Nblock+1)] for j in range(Nconf)]
    esln_ave = [ 0 for i in range(Nblock+1)]
    for l in range(Nst, Ned+1):
        for k in range(1, Nblock+1):
            inpf_sln = inpdir + 'ERmod_{0:04d}/soln/engsln.{1:02d}'.format(l,k)
            print(inpf_sln)
            with open(inpf_sln) as f:
                esln[l-1][k-1] = [line.split() for line in f]
#                print(esln[l-1][k-1])
                esln_ave[k-1] = ['' for line in range(len(ene_coord))]

#    print(esln[0])
    sum_ = [[0 for line in range(num_lines)] for line in range(Nblock)]
    ave = [[0 for line in range(num_lines)] for line in range(Nblock)]
    print('==========\n\n==========')
    for ie in range(num_lines):
        for k in range(1, Nblock+1):
            sum_[k-1][ie] = 0.0e0
            for l in range(Nst, Ned):
                print('l:', l, 'k:', k)
                print(esln[l-1][k-1][ie][2])
                sum_[k-1][ie] += float(esln[l-1][k-1][ie][2]) 
                ave[k-1][ie] = sum_[k-1][ie]/ Nconf
                
            esln_ave[k-1][ie] = [ene_coord[ie], '1', ave[k-1][ie]]
            print(esln_ave[k-1][ie])
    for k in range(1, Nblock+1):
        outf = inpdir + 'aveERmod_{0:02d}/soln/engsln.{1:02d}'.format(1,k)
        shutil.move(outf, outf+'.bak')
        with open(outf, 'a+') as f:
            for ie in range(num_lines):
                f.write(esln_ave[k-1][ie][0]+'\t'+esln_ave[k-1][ie][1]+'\t'+'{0:12.9e}'.format(esln_ave[k-1][ie][2])+'\n')
