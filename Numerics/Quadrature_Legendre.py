#!/usr/bin/env python3
import numpy
import math
import scipy.special as spsp
import scipy.linalg as spla
import numpy.linalg as npla

if (__name__ == '__main__'):

    N = 4

    PN=spsp.legendre(N)
    print(PN.roots)
    print(type(PN.roots))
    print('ABSCISSAS',numpy.sort(PN.roots))

    X = numpy.sort(PN.roots)
    M = numpy.zeros([N,N],dtype=numpy.float64)
    for i in range(N):
        P = spsp.legendre(i)
        p = numpy.zeros([N])
        for j in range(N):
            M[i,j] = P(X[j]) #/(math.sqrt(2/(2*i+1)))

    M2=numpy.matmul(M,M.transpose())
    #W,vr=spla.eig(M)
    print(M2)
    W2,vr=spla.eigh(M2)
    print()
    print('EIGVALS',W2)
    print()
    print('EIGVECS')
    print(vr)
    print()
    print('ORTHOGONALITY CHECK')
    print(numpy.matmul(vr,vr.transpose()))


