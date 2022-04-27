import numpy as np
from scipy.sparse import csr_matrix
import multiprocessing as mp


def get_index(ind):
    return np.where(p_address[:] == ind)[0]


def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)),
                      shape=(data.max() + 1, data.size))

def get_indices_sparse(data):
    M = compute_M(data)
    return [np.unravel_index(row.data, data.shape) for row in M]


rN = 2
n_mesh = rN**3
N = 32
L = 2.0
dx = 1.0
meshfilter = np.array([1,rN,rN**2])
data = np.arange(n_mesh)
print(data)

pos = np.random.rand(N,3)*L

# p_address
p_address = np.sum((pos/dx).astype(int)*meshfilter,axis=1)
print(p_address)

# multi process:
p = mp.Pool()
address_list = p.map(get_index, data)
print('address 0:')
print(address_list)

# multi process:
n_process = int(N/n_mesh)
p = mp.Pool(processes=n_process)
pair_list = p.map(get_indices_sparse, np.array([p_address]))
print('address 1:')
print(pair_list)
