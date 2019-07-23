import ctypes
import numpy as np

def real_fortran(i, j, k, box, rcut, kappa, charge, A, e_real):
    f = np.ctypeslib.load_library("librealfort.so", ".")
    f.ureal_.argtypes = [
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
    f.ureal_.restype = ctypes.c_void_p

    fi = ctypes.byref(ctypes.c_int32(i))
    fj = ctypes.byref(ctypes.c_int32(j))
    fk = ctypes.byref(ctypes.c_int32(k))
    box = ctypes.byref(ctypes.c_double(box))
    rcut = ctypes.byref(ctypes.c_double(rcut))
    kappa = ctypes.byref(ctypes.c_double(kappa))

    f.ureal_(fk, fj, fi, box, rcut, kappa, charge, A, e_real)

i =  5
a =  4
v =  3
box = 3.0
rcut = 1.20
kappa = 1.0 

charge = np.zeros((i, a))
e_real = np.zeros((i,i))
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

real_fortran(i, a, v, box, rcut, kappa, charge, A, e_real)


print(e_real)
