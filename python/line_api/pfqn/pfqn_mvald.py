from numpy import *
from ..util import *
def pfqn_mvald(L = None,N = None,Z = None,mu = None,stabilize = None):
    # [XN,QN,UN,CN,LGN]=PFQN_MVALD(L,N,Z,MU)

    # stabilize ensures that probabilities do not become negative
    if len(varargin) < 5:
        stabilize = True

    warn = True
    isNumStable = True
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]

    Xs = zeros((R,prod(N + 1)))

    pi = ones((M,sum(N) + 1,prod(N + 1)))

    WN = zeros((M,R))
    n = pprod(N)

    lGN = 0
    while n != - 1:

        WN = 0 * WN
        for s in arange(0,R+1):
            if N[s] > 0:
                for i in arange(0,M+1):
                    WN[i,s] = 0
                    for k in arange(0,sum(n)+1):
                        WN[i,s] = WN(i,s) + (L(i,s) / mu(i,k)) * k * pi(i,(k - 1) + 1,hashpop(oner(n,s),N))
                Xs[s,hashpop[n,N]] = N[s] / (Z[s] + sum(WN(:,s)))
            # compute pi(k|n)
            for k in arange(0,sum(n)+1):
                for i in arange(0,M+1):
                    pi[i,[k] + 1,hashpop[n,N]] = 0
                for s in arange(0,R+1):
                    if N[s] > 0:
                        for i in arange(0,M+1):
                            pi[i,[k] + 1,hashpop[n,N]] = pi(i,[k] + 1,hashpop(n,N)) + (L(i,s) / mu(i,k)) * Xs(s,hashpop(n,N)) * pi(i,(k - 1) + 1,hashpop(oner(n,s),N))
            # compute pi(0|n)
            for i in arange(0,M+1):
                p0 = 1 - sum(pi(i,(arange(0,sum(n)+1)) + 1,hashpop(n,N)))
                if p0 < 0:
                    if warn:
                        line_warning(mfilename,'MVA-LD is numerically unstable on this model, forcing all probabilities to be non-negative.')
                        #                N
                        warn = False
                        isNumStable = False
                    if stabilize:
                        pi[i,[0] + 1,hashpop[n,N]] = eps
                    else:
                        pi[i,[0] + 1,hashpop[n,N]] = p0
                else:
                    pi[i,[0] + 1,hashpop[n,N]] = p0
            last_nnz = find(n > 0,1,'last')
            if sum(n(arange(0,last_nnz - 1+1))) == sum(N(arange(0,last_nnz - 1+1))) and sum(n(arange((last_nnz + 1),R+1))) == 0:
                logX = log(Xs(last_nnz,hashpop(n,N)))
                #hashpop(n,N)
                if not len(logX)==0 :
                    lGN[end() + 1] = lGN(end()) - logX
            n = pprod(n,N)


        X = Xs(:,hashpop(N,N))
        XN = transpose(X)
        pi = pi(:,:,hashpop(N,N))
        QN = multiply(WN,matlib.repmat(XN,M,1))
        #UN = repmat(XN,M,1) .* L;
        UN = 1 - pi(:,1)
        CN = N / XN - Z

        return XN,QN,UN,CN,lGN,isNumStable,pi