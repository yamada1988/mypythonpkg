import numpy as np
import cupy as cp

V = 3.048**3.0E0 #(nm^3)
kBT = 2.479 #(kJ/mol)
alpha = kBT/V # (kJ/mol/nm^3)
alpha = alpha * 1.602*10**(6.0E0) #(Pa)
alpha = alpha * 10**(-5.0E0) #(bar)


fname = 'sigma.xvg'

with open(fname, 'rt') as f:
    sigma = np.array([float(line.split()[1]) for line in f if not line.startswith('#')])
    f.seek(0)
    t = np.array([float(line.split()[0]) for line in f if not line.startswith('#')])

N = len(t)
N = 100
G = cp.zeros(N, dtype=np.float32)
Z = cp.zeros(N)
sigma = cp.asarray(sigma[:N], dtype=np.float32)
it0 = 0
data = cp.array(range(N), dtype=np.int32)

g = cp.zeros(N, dtype=np.float32)
mat_add_kernel = cp.ElementwiseKernel(
            in_params='int32 data, raw float32 P, int32 N',
            out_params='raw float32 g',
            operation=\
            '''
            for (int it0 = 0; it0 < N-i; it0++){
              g[i] += P[i+it0]*P[it0];
            }
            g[i] /= float(N-i)
            ''',
            name='mat_add_kernel')

mat_add_kernel(data, sigma, N, g)
g = cp.asnumpy(g*10**5.0/alpha)
for i in range(N):
    print(t[i], g[i])
