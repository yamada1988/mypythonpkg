import ctypes
import numpy as np

def calc_fortran(i, j, k, A):
    f = np.ctypeslib.load_library("libfort.so", ".")
    f.calc_.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64)
        ]
    f.calc_.restype = ctypes.c_void_p

    fi = ctypes.byref(ctypes.c_int32(i))
    fj = ctypes.byref(ctypes.c_int32(j))
    fk = ctypes.byref(ctypes.c_int32(k))

    f.calc_(fi, fj, fk, A)



i =  5
a =  4
v =  3
box = 3.0
rcut = 1.20

A = box*np.random.rand(i, a, v)

for ii in range(i):
        for ia in range(a):
                print(ii, ia, A[ii,ia,:])

calc_fortran(i, a, v, A)
