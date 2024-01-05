from numpy import *
import numpy.linalg as linalg

def pfqn_kt(L = None,N = None,Z = None):
    # Knessl-Tier asymptotic expansion

    if len(L)==0 or len(N)==0 or sum(N) == 0:
        G = 1
        return G,lG,X,Q

    if Z is None:
        Z = N * 0

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    Ntot = sum(N)
    beta = ones((1,R))
    for r in arange(0,R+1):
        beta[r] = N[r] / Ntot

    X,Q = pfqn_aql(L,N,Z)
    delta = eye(R,R)
    for i in arange(0,R+1):
        for j in arange(0,R+1):
            SK = 0
            for k in arange(0,M+1):
                SK = SK + X[i] * X[j] * L[k,i] * L[k,j] / (1 - sum(transpose(X) * transpose(L[k,:]))) ** 2
            C[i,j] = delta[i,j] * beta[i] + (1 / Ntot) * SK

    Den = 1
    for k in arange(0,M+1):
        Den = Den * amax(1e-06,(1 - sum(transpose(X) * transpose(L[k,:]))))

        #G=(2*pi).^(-R/2)/sqrt(Ntot^R*det(C))*exp(-Ntot*beta*log(X)')/Den;
        lG = log((2 * pi) ** (- R / 2) / sqrt(Ntot ** R * linalg.det(C))) + (- Ntot * beta * transpose(log(X))) - log(Den)
        G = exp(lG)
        return G,lG,X,Q
