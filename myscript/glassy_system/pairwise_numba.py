# coding: utf-8
import numpy
from numba import double
from numba.decorators import jit
import time
import traceback

def pairwise(X, Y, D):
    M0 = X.shape[0]
    M1 = Y.shape[0]
    N =  X.shape[1]
    for i in range(M0):
        for j in range(M1):
            d = 0.0
            for k in range(N):
                tmp = X[i, k] - Y[j, k]
                tmp -= int(tmp)*tmp
                d += tmp * tmp
            D[i, j] = numpy.sqrt(d)

    return D


@jit
def pairwise_numba(X, Y, D):
    M0 = X.shape[0]
    M1 = Y.shape[0]
    N =  X.shape[1]
    for i in range(M0):
        for j in range(M1):
            d = 0.0
            for k in range(N):
                tmp = X[i, k] - Y[j, k]
                tmp -= int(tmp)*tmp
                d += tmp * tmp
            D[i, j] = numpy.sqrt(d)

    return D

@jit('f8[:, :](f8[:, :], f8[:,:], f8[:, :])', nopython=True)
def pairwise_numba2(X, Y, D):
    M0 = X.shape[0]
    M1 = Y.shape[0]
    N = X.shape[1]
    for i in range(M0):
        for j in range(M1):
            d = 0.0
            for k in range(N):
                tmp = X[i, k] - Y[j, k]
                tmp -= int(tmp)*tmp
                d += tmp * tmp
            D[i, j] = numpy.sqrt(d)

    return D

if __name__ == '__main__':

    N = 120
    t = time.time()
    for it in range(16**3):
        X = numpy.random.random((N, 3))
        Y = numpy.random.random((8,3))
        D = numpy.empty((N, 8))
        pairwise(X, Y, D)
    print "vanilla:", time.time() - t

    t = time.time()
    for it in range(16**3):
        X = numpy.random.random((N, 3))
        Y = numpy.random.random((8,3))
        D = numpy.empty((N, 8))
        pairwise_numba(X, Y, D)
    print "numba:", time.time() - t

    t = time.time()
    for it in range(16**3):
        X = numpy.random.random((N, 3))
        Y = numpy.random.random((8,3))
        pairwise_numba2(X, Y, D)
    print "numba2:", time.time() - t
