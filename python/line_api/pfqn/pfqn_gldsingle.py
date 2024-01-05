from numpy import *

def pfqn_gldsingle(L = None,N = None,mu = None,options = None):
    # G=PFQN_GLDSINGLE(L,N,MU)

    if len(varargin) < 4:
        options = []

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    assert R > 1, 'Multiclass model detected. pfqn_gldsingle is for single class models.'

    g = L(1,1) * 0
    for n in arange(0,N+1):
        g[0 + 1,n + 1,1 + 1] = 0

    for m in arange(0,M+1):
        for tm in arange(0,(N + 1)+1):
            g[m + 1,0 + 1,tm + 1] = 1
        for n in arange(0,N+1):
            for tm in arange(0,(N - n + 1)+1):
                g[m + 1,n + 1,tm + 1] = g(m - 1 + 1,n + 1,1 + 1) + L(m) * g(m + 1,n - 1 + 1,tm + 1 + 1) / mu(m,tm)

    G = g(M + 1,N + 1,1 + 1)
    lG = log(G)
    return lG,G