import ctypes
import numpy as np

def lj_fortran(i, j, k, box, rcut, sigma, eps, A, e_lj):
    f = np.ctypeslib.load_library("libljfort.so", ".")
    f.lj_.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)
        ]
    f.lj_.restype = ctypes.c_void_p

    fi = ctypes.byref(ctypes.c_int32(i))
    fj = ctypes.byref(ctypes.c_int32(j))
    fk = ctypes.byref(ctypes.c_int32(k))
    box = ctypes.byref(ctypes.c_double(box))
    rcut = ctypes.byref(ctypes.c_double(rcut))

    f.lj_(fi, fj, fk, box, rcut, sigma, eps, A, e_lj)






i =  1000
a =  4
v =  3
box = 3.0
rcut = 1.20

sigma = np.zeros((i, a))
eps = np.zeros((i, a))
e_lj = np.zeros((i,i))
for im in range(i):
    for ia in range(a):
        sigma[im,ia] = 0.40*(ia+1)
        eps[im, ia] =  0.40*(ia+1)

A = box*np.random.rand(i, a, v)

#for ii in range(i):
#    for ij in range(a):
#        for ji in range(ii+1,i):
#            for jj in range(a):
#                r = A[ii,ij] - A[ji,jj]
#                r -= box*np.round(r/box)
#                print(ii, ji, ij, jj, r)

import time
t0 = time.time()
lj_fortran(i, a, v, box, rcut, sigma, eps, A, e_lj)
print('fortran:', time.time()-t0)

print(e_lj)
