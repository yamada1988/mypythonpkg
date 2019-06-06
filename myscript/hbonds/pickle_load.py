import pickle
from scipy.sparse import lil_matrix, csr_matrix

tN = 700
sij = []
for it in range(tN):
    with open('sij/sij_{0:05d}.pickle'.format(it), 'rb') as f:
        s = pickle.load(f)
    sij.append(s)

cij = [0.0e0 for i in range(tN)]
for i0, s0 in enumerate(sij):
    print(i0, s0.sum())
    s0_csr = s0.tocsr()
    for it, s in enumerate(sij[i0:]):
        s_csr = s.tocsr()
        temp = s0_csr.multiply(s_csr)
        cij[it] += temp.sum()

for it in range(tN):
    cij[it] /= (tN - it)

sum_ave = [0.0e0 for i in range(tN)]
for it, s0 in enumerate(sij):
    s0_csr = s0.tocsr()
    temp = s0_csr.multiply(s0_csr)
    sum_ave[it] = temp.sum()


dt = 0.20e0 #ps
with open('chb.dat', 'wt') as f:
    for it in range(tN):
        cij[it] /= (sum(sum_ave[:tN-it])/(tN-it))
        print(it, cij[it])
        f.write('{0:5.3f}\t{1:11.10f}\n'.format(it*dt, cij[it]))
