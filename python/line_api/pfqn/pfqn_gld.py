from numpy import *

from ..api.pfqn.pfqn_gldsingle import pfqn_gldsingle
from ..util import *

def pfqn_gld(L = None,N = None,mu = None,options = None):
    # [G,LG]=PFQN_GLD(L,N,MU,OPTIONS)

    # G=pfqn_gld(L,N,mu)
    # mu: MxN matrix of load-dependent rates
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    lambda_ = zeros((1,R))
    if M == 1:
        lG = factln(sum(N)) - sum(factln(N)) + N(L > 0) * transpose(log(L(L > 0))) - sum(log(mu(1,arange(0,sum(N)+1))))
        G = exp(lG)
        return G,lG

    if R == 1:
        lG,G = pfqn_gldsingle(L,N,mu)
        return G,lG

    if len(L)==0:
        G = 0
        lG = - Inf
        return G,lG

    if mu is None:
        mu = ones((M,sum(N)))

    if options is None:
        options = SolverNC.defaultOptions

    isLoadDep = False
    isInfServer = []
    for i in arange(0,M+1):
        if amin(mu[i,arange(0,sum(N)+1)]) == logical_and(1,amax(mu[i,arange(0,sum(N)+1)])) == 1:
            isInfServer[i] = False
            continue
        else:
            if all(arange(mu[i,arange(0,sum(N)+1)]== 1,sum(N)+1)):
                isInfServer[i] = True
                continue
            else:
                isInfServer[i] = False
                isLoadDep = True

    if not isLoadDep :
        # if load-independent model then use faster pfqn_gmva solver
        Lli = L[nonzero(not isInfServer ),:]
        if len(Lli)==0:
            Lli = 0 * N
        Zli = L[nonzero(isInfServer),:]
        if len(Zli)==0:
            Zli = 0 * N
        options.method = 'exact'
        lG = pfqn_nc(lambda_,Lli,N,sum(Zli, 0),options)
        G = exp(lG)
        return G,lG

    G = 0
    if M == 0:
        G = 0
        lG = log(G)
        return G,lG

    if sum(N == zeros((1,R))) == R:
        G = 1
        lG = log(G)
        return G,lG

    if R == 1:
        G = pfqn_gldsingle(L,N,mu)
        lG = log(G)
        return G,lG

    G = G + pfqn_gld(L[arange(0,(M - 1)+1),:],N,mu[arange(0,(M - 1)+1),:], options)
    for r in arange(0,R+1):
        if N[r] > 0:
            if R > 1:
                N_1 = oner(N,r)
            else:
                N_1 = N - 1
            G = G + (L[M,r] / mu[M,1]) * pfqn_gld(L,N_1,pfqn_mushift(mu,M))

    lG = log(G)
    return G,lG
