import numpy as np
import sys

args = sys.argv
dirname = args[1]

histss = []
histsr = []
blocknum = 10
for ir in range(1,blocknum+1):
    rfname = dirname + '/aveERmod_confs/block{0:02d}/histogram_flcuv_refs.tt'.format(ir)
    sfname = dirname + '/aveERmod_confs/block{0:02d}/histogram_flcuv_soln.tt'.format(ir)
    des = np.loadtxt(sfname)[:,0]
    histss.append(np.loadtxt(sfname)[:,1])
    der = np.loadtxt(rfname)[:,0]
    histsr.append(np.loadtxt(rfname)[:,1])
    N = len(der)

histss_ave = np.zeros(N)
histsr_ave = np.zeros(N)
histss_err = np.zeros(N)
histsr_err = np.zeros(N)

de1 = der[1]-der[0]
de2 = 5.0
de3 = 100000000.0

#print(len(histss[0]))
for i in range(N): 
    hs_s = []
    hs_r = []
    for ir in range(blocknum):
        hs_s.append(histss[ir][i])
        hs_r.append(histsr[ir][i])
    hs_s = np.array(hs_s)
    hs_r = np.array(hs_r)
    #print(hs_s,hs_r)
    histss_ave[i] = np.mean(hs_s)
    histsr_ave[i] = np.mean(hs_r)
    #print(histss_ave[i], histsr_ave[i])
    histss_err[i] = np.std(hs_s)
    histsr_err[i] = np.std(hs_r)

    #print(histss_err[i], histsr_err[i])


outf1 = dirname + '/histogram_flcuv_soln_error.tt'
with open(outf1, 'wt') as f:
     l = '#e\thist\terr\n'
     f.write(l)
     for i in range(N):
         if i%10 == 0:
             l = "{0:7.5e}\t{1:9.7e}\t{2:9.7e}\n".format(des[i], histss_ave[i], histss_err[i])
             f.write(l)

outf2 = dirname + '/histogram_flcuv_refs_error.tt'
with open(outf2, 'wt') as f:
     l = '#e\thist\terr\n'
     f.write(l)
     for i in range(N):
         if i%10 == 0:
             l = "{0:7.5e}\t{1:9.7e}\t{2:9.7e}\n".format(der[i], histsr_ave[i], histsr_err[i])
             f.write(l)
