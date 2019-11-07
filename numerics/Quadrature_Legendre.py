#!/usr/bin/env python3
import numpy
import scipy.special as spsp
import scipy.linalg as spla
import numpy.linalg as npla

if (__name__ == '__main__'):

    N = 5

    PN=spsp.legendre(N)
    print(PN.roots)
    print(type(PN.roots))
    print(numpy.sort(PN.roots))

    X = numpy.sort(PN.roots)
    M = numpy.zeros([N,N],dtype=numpy.float64)
    Msq = numpy.zeros([N,N],dtype=numpy.float64)
    for i in range(N):
        P = spsp.legendre(i)
        p = numpy.zeros([N])
        for j in range(N):
            M[i,j] = P(X[j])
            Msq[i,j] = (P(X[j]))**2

    M2=numpy.matmul(M,M.transpose())
    W,vr=spla.eig(M)
    print(W)
    W2,vr=spla.eig(M2)
    print(W2)

    Ones = numpy.ones(N)
    w = npla.solve(Msq,Ones)
    print(w)


