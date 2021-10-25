#!/usr/bin/env python3
'''
Implementation of Numerov algorithm

D. A. Brue

'''
import math
import numpy as np
import logging


def plotSolution(X,V,Fs,Fl,Fr,T):
    import matplotlib.pyplot as plt

    fig0 = plt.figure()
    ax0 = fig0.add_subplot(1,1,1)


    plt.plot(X,V,'k-')

    if (type(Fl) == list ):
        for i in range(len(Fs)):
            plt.plot(X,Fs[i],'k-',label='Solution')
        #for i in range(len(Fl)):
        #    plt.plot(X,Fl[i],'b-',label='Left Sol.')
        #for i in range(len(Fr)):
        #    plt.plot(X,Fr[i],'r-',label='Right Sol.')
    else:
        plt.plot(X,Fl,'b-',label='Left Sol.')
        plt.plot(X,Fr,'r-',label='Right Sol.')
    plt.ylim([-math.pi,3*math.pi])
    plt.show()

def matchLR(X,Fl,Fr,matchi = None):
    '''
    tol: optional argument for tolerance in matching algorithm
    matchi: int: index of X to perform the matching. If None, midpoint
    '''

    if (not matchi):
        # if a matching index point is not supplied, try midpoint
        matchi = int((len(X)/2))
    elif (matchi < 1 or matchi >= len(X)-1):
        print('Bad matching index, out of bounds')
        exit()

    # Match value, check derivative
    tx = np.array([X[matchi-1],X[matchi],X[matchi+1]])
    tl = np.array([Fl[matchi-1],Fl[matchi],Fl[matchi+1]])
    tr = np.array([Fr[matchi-1],Fr[matchi],Fr[matchi+1]])
    ratio = tl[1]/tr[1]
    tr *= ratio
    ml = tl[2]-tl[0]
    mr = tr[2]-tr[0]
    mdiff = ml-mr
    mrat = abs(ml/mr) - 1
    #mrat = (ml/mr) - 1
    mdet = (tr[0]*tl[2]-tr[2]*tl[0])

    return mdiff, mrat, mdet
    

def int1dLR(X, V, E, matchi=None):
    '''

    Integrate in 1D for a given function V such that the solution function
    gives 

        (d/dx)^2 F(x) + V(X)*F(x) = E*F(x)

    X: np.array: 1d spatial array
    V: np.array: potential curve
    E: float: system energy, potential eigenvalue

    Output:
    F(x): np.array, X-like

    See accompanying documentation for algorithm derivation
    '''

    log = logging.getLogger(__name__)

    initVal = 1.0e-3
    h = X[1]-X[0]
    N = len(X)

    # Build T array
    T = np.zeros_like(X)
    T = E - V

    if (not matchi):
        matchi = int(N/2)

    if (T[0] > 0.0 or T[-1] > 0.0):
        log.warning('WARNING: Numerov 1D, end-points do not lead to local solution')
        print(T[0],T[-1],E,V[0],V[-1])
        exit()

    classical_bound_index_L = 0
    for i in range(N):
        if (T[i] > 0.0):
            classical_bound_index_L = i
            break
    classical_bound_index_R = 0
    for i in range(N-1,-1,-1):
        if (T[i] > 0.0):
            classical_bound_index_R = i
            break

    Fl = np.zeros_like(X)
    Fr = np.zeros_like(X)

    # Left-to-Right propagation
    Fl[1] = initVal
    for i in range(2,classical_bound_index_R):
        Fl[i] = ( Fl[i-1]*(2+h**2*T[i-1]*5/6)- Fl[i-2]*(1-h**2*T[i-2]/12))/(1-T[i]*h**2/12)
        if (abs(Fl[i]) > 1.0e9):
            Fl *= 1/1.0e8

    # Right-to-Left propagation
    Fr[-2] = initVal
    for i in range(N-3,classical_bound_index_L,-1):
        Fr[i] = ( Fr[i+1]*(2+h**2*T[i+1]*5/6)- Fr[i+2]*(1-h**2*T[i+2]/12))/(1-T[i]*h**2/12)
        if (abs(Fr[i]) > 1.0e9):
            Fr *= 1/1.0e8

    F = np.zeros_like(X)
    F[:matchi] = Fl[:matchi]
    F[matchi:] = Fr[matchi:]
    norm = math.sqrt(np.sum(F**2) * (X[1]-X[0]))
    if (norm > 0):
        F *= 1/norm

    return F,Fl,Fr,T


####################################################################################################
def TEST0():
    print('Running Numerov TEST0: Components')
    import matplotlib.pyplot as plt

    # Solve harmonic potential, V = (1/2)*k*x^2

    nPts = 1001

    matchi = int(nPts/2)
    
    k = 10.0

    xMin = -10*math.pi
    xMax =  10*math.pi

    nE = 4901
    Emin = 1.0
    Emax = 50.0
    Ens = np.linspace(Emin,Emax,nE)

    X = np.linspace(xMin,xMax,nPts)
    V = k*X**2/2

    Fs = []
    Fls = []
    Frs = []
    dets = []
    rats = []
    difs = []
    for i in range(nE):
        E = Ens[i]
        print('Energy ',E)
        F,Fl,Fr,T = int1dLR(X,V,E)
        Fs.append(F)
        Fls.append(Fl)
        Frs.append(Fr)
        mdif,mrat,mdet = matchLR(X,Fl,Fr,matchi=matchi)
        dets.append(mdet)
        difs.append(mdif)
        rats.append(mrat)

    # Clip large values
    if (True):
        clip = 1.0e2
        for i in range(len(rats)):
            if (rats[i] > clip):
                rats[i] = clip
            if (rats[i] < -clip):
                rats[i] = -clip
        for i in range(len(difs)):
            if (difs[i] > clip):
                difs[i] = clip
            if (difs[i] < -clip):
                difs[i] = -clip
        for i in range(len(dets)):
            if (dets[i] > clip):
                dets[i] = clip
            if (dets[i] < -clip):
                dets[i] = -clip

    #rtn = plotSolution(X,V,Fs,Fls,Frs,T)

    fig0=plt.figure()
    ax0 = fig0.add_subplot(1,1,1)
    plt.plot(Ens,difs,'r*',label='Differences')
    plt.plot(Ens,rats,'g.',label='Ratios')
    plt.plot(Ens,dets,'b.',label='Determinant')
    #plt.ylim([-10,10])
    plt.show()

    return 0

####################################################################################################
def TEST1():
    print('Running Numerov TEST1: Multiproc')
    import matplotlib.pyplot as plt
    import multiprocessing as mp

    # Solve harmonic potential, V = (1/2)*k*x^2

    # INIT VALS -------------------------------
    nPts = 1001
    matchi = int(nPts/2)
    k = 10.0
    xMin = -10*math.pi
    xMax =  10*math.pi
    nE = 4901
    Emin = 1.0
    Emax = 50.0
    Ens = np.linspace(Emin,Emax,nE)
    X = np.linspace(xMin,xMax,nPts)
    V = k*X**2/2

    # Create the mulitproc pool, typically number of cpus
    funcMap = mp.Pool()
    # Share the potential and space arrays between processes (removes massive duplicates)
    Vshare = mp.Array(V)
    Xshare = mp.Array(X)

    def runNum1D(X,V,E):
        F,Fl,Fr,T = int1dLR(X,V,E)
        mdif,mrat,mdet = matchLR(X,Fl,Fr,matchi=matchi)
        return [E,F,mdif,mrat,mdet]

    Fs = []
    Fls = []
    Frs = []
    dets = []
    rats = []
    difs = []
    for i in range(nE):
        E = Ens[i]
        print('Energy ',E)
        F,Fl,Fr,T = int1dLR(X,V,E)
        Fs.append(F)
        Fls.append(Fl)
        Frs.append(Fr)
        mdif,mrat,mdet = matchLR(X,Fl,Fr,matchi=matchi)
        dets.append(mdet)
        difs.append(mdif)
        rats.append(mrat)

    # Clip large values
    if (True):
        clip = 1.0e2
        for i in range(len(rats)):
            if (rats[i] > clip):
                rats[i] = clip
            if (rats[i] < -clip):
                rats[i] = -clip
        for i in range(len(difs)):
            if (difs[i] > clip):
                difs[i] = clip
            if (difs[i] < -clip):
                difs[i] = -clip
        for i in range(len(dets)):
            if (dets[i] > clip):
                dets[i] = clip
            if (dets[i] < -clip):
                dets[i] = -clip

    #rtn = plotSolution(X,V,Fs,Fls,Frs,T)

    fig0=plt.figure()
    ax0 = fig0.add_subplot(1,1,1)
    plt.plot(Ens,difs,'r*',label='Differences')
    plt.plot(Ens,rats,'g.',label='Ratios')
    plt.plot(Ens,dets,'b.',label='Determinant')
    #plt.ylim([-10,10])
    plt.show()

    return 0
if (__name__ == '__main__'):
    ''' Run some test cases '''
    #rtn=TEST0()
    rtn=TEST1()
