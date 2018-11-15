import cupy as cp
import numpy as np
import time

with open('benchmark.dat', 'w') as f:
    f.write('# totalatm\tcupy\tnumpy\n')
rN = 64
for N in [5000, 10000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000]:
    print('============\nN={0:9d}\n============'.format(N))
    print('cupy:')
    r = cp.asarray(np.random.randint(0,rN,(N,3), dtype=np.int32))
    n = cp.asnumpy(r)
    w = cp.asarray(np.random.rand(N), dtype=np.float32)
    m = cp.asnumpy(w)
    y = cp.zeros((rN,rN,rN), dtype=cp.float32)
    g = cp.asnumpy(y)
    t0 = time.time()
    r0 = r[:,0]
    r1 = r[:,1]
    r2 = r[:,2]
    cp.ElementwiseKernel(
            'int32 r0, int32 r1, int32 r2, float32 w',
            'raw float32 y',
            '''
            int ind[] = {r0,r1,r2};
            atomicAdd(&y[ind], w);
            '''
            )(r0,r1,r2,w,y)

    dt_cp = time.time()-t0
    print('time:', dt_cp)
    y = cp.asnumpy(y)

    print('numpy:')
    t1 = time.time()
    m, ax = np.histogramdd(n, bins=(rN, rN, rN), range=((0, rN), (0, rN), (0, rN)), weights=m)

    dt_np = time.time() - t1
    print('time:', dt_np) 
    with open('benchmark.dat','a+') as f:
        f.write('{0:09d}\t{1:7.5f}\t{2:7.5f}\n'.format(N, dt_cp, dt_np))

