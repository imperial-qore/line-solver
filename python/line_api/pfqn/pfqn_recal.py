from numpy import *

def pfqn_recal(L = None,N = None,Z = None,m0 = None):
    # [G,logG]=PFQN_RECAL(L,N,Z,M0)
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    Ntot = sum(N)
    G_1 = ones((1,nchoosek(Ntot + M - 1,Ntot)))
    G = G_1
    if len(varargin) < 4:
        m0 = ones((1,M))

    if len(varargin) < 3 or sum(Z) == 0:
        I_1 = multichoose(M,Ntot)
        n = 0
        for r in arange(0,R+1):
            for nr in arange(0,N[r]+1):
                n = n + 1
                I = multichoose(M,(Ntot + 1) - (n + 1))
                for i in arange(0,I.shape[0]+1):
                    m = I[i,:]
                    G[i] = (0)
                    for j in arange(0,M+1):
                        m[j] = m[j] + 1
                        G[i] = G[i] + (m[j] + m0[j] - 1) * L[j,r] * G_1(matchrow(I_1,m)) / nr
                        m[j] = m[j] - 1
                I_1 = I
                G_1 = G

    G = G(1)
    lG = log(G)
    return G,lG