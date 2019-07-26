import ctypes
import numpy as np

def corr_fortran(i, j, k, box, rcut, kappa, charge, A, e_corr):
    f = np.ctypeslib.load_library("libcorrfort.so", ".")
    f.ucorr_.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)
        ]
    f.ucorr_.restype = ctypes.c_void_p

    fi = ctypes.byref(ctypes.c_int32(i))
    fj = ctypes.byref(ctypes.c_int32(j))
    fk = ctypes.byref(ctypes.c_int32(k))
    box = ctypes.byref(ctypes.c_double(box))
    rcut = ctypes.byref(ctypes.c_double(rcut))
    kappa = ctypes.byref(ctypes.c_double(kappa))

    f.ucorr_(fi, fj, fk, box, rcut, kappa, charge, A, e_corr)

i =  5
a =  4
v =  3
box = 3.0
rcut = 1.20
kappa = 1.0 

charge = np.zeros((i, a))
e_corr = np.zeros((i,i))
for im in range(i):
    for ia in range(a):
        charge[im,ia] = 0.40*(ia+1)

A = box*np.random.rand(i, a, v)


for ii in range(i):
    for ij in range(a):
        for ji in range(ii+1,i):
            for jj in range(a):
                r = A[ii,ij] - A[ji,jj]
                r -= box*np.round(r/box)
                print(ii, ji, ij, jj, r)

corr_fortran(i, a, v, box, rcut, kappa, charge, A, e_corr)


print(e_corr)
