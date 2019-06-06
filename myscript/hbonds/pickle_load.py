import pickle
from scipy.sparse import lil_matrix, csr_matrix


tqdmflag = 0
try:
    from tqdm import trange, tqdm
except:
    tqdmflag = 1

tqdmflag = 1
tN = 2000
dt = 0.20e0 #ps
sij = []
if tqdmflag == 0:
    for it in trange(tN, desc = 'load pickled files'):
        with open('sij/sij_{0:05d}.pickle'.format(it), 'rb') as f:
            s = pickle.load(f)
        sij.append(s)
elif tqdmflag == 1:
    for it in range(tN):
        with open('sij/sij_{0:05d}.pickle'.format(it), 'rb') as f:
            s = pickle.load(f)
        sij.append(s)


cij = [0.0e0 for i in range(tN)]
if tqdmflag == 0:
    for i0 in trange(tN, desc='hbond lifetime calculation'):
        s0 = sij[i0]
        s0_csr = s0.tocsr()
        for it, s in enumerate(sij[i0:]):
            s_csr = s.tocsr()
            temp = s0_csr.multiply(s_csr)
            cij[it] += temp.sum()/s0_csr.sum()
else:
   for i0, s0 in enumerate(sij):
       s0_csr = s0.tocsr()
       for it, s in enumerate(sij[i0:]):
           s_csr = s.tocsr()
           temp = s0_csr.multiply(s_csr)
           cij[it] += temp.sum()/s0_csr.sum()


for it in range(tN):
    cij[it] /= (tN - it)

#sum_ave = [0.0e0 for i in range(tN)]
#for it, s0 in enumerate(sij):
#    s0_csr = s0.tocsr()
#    temp = s0_csr.multiply(s0_csr)
#    sum_ave[it] = temp.sum()


with open('CHB.dat', 'wt') as f:
    for it in range(tN):
        cij[it] /= cij[0]
        print(it, cij[it])
        f.write('{0:5.3f}\t{1:11.10f}\n'.format(it*dt, cij[it]))
