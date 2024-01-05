from numpy import *
from ..util import *
def pfqn_panacea(L = None,K = None,Z = None):
    # [GN,LGN]=PFQN_PANACEA(L,K,Z)

    # K = population vector
    q,p = L.shape
    if len(varargin) == 2 or len(Z)==0:
        Z = K * 0 + 1e-08

    if logical_or(len(L)==0,sum(L, 0)) == zeros((1,p)):
        lGn = - sum(factln(K)) + sum(multiply(K,log(sum(Z, 0))))
        Gn = exp(lGn)
        return Gn,lGn

    r = L / matlib.repmat(Z,q,1)
    N = amax(1.0 / r)
    beta = K / N
    gamma = r * N
    alpha = 1 - K * transpose[r]
    gammatilde = gamma / matlib.repmat(transpose(alpha),1,p)
    if amin(alpha) < 0:
        #    line_warning(mfilename,'Model is not in normal usage');
        Gn = NaN
        lGn = NaN
        return Gn,lGn

    A0 = 1
    A1 = 0
    for j in arange(0,p+1):
        m = zeros((1,p))
        m[j] = 2
        A1 = A1 - beta[j] * pfqn_ca(gammatilde,m)

    A2 = 0
    for j in arange(0,p+1):
        m = zeros((1,p))
        m[j] = 3
        A2 = A2 + 2 * beta[j] * pfqn_ca(gammatilde,m)
        m = zeros((1,p))
        m[j] = 4
        A2 = A2 + 3 * beta[j] ** 2 * pfqn_ca(gammatilde,m)
        for k in setdiff1d(arange(0,p+1),j):
            m = zeros((1,p))
            m[j] = 2
            m[k] = 2
            A2 = A2 + 0.5 * beta[j] * beta[k] * pfqn_ca(gammatilde,m)

    if 0:
        A3 = 0
        for j in arange(0,p+1):
            m = zeros((1,p))
            m[j] = 4
            A3 = A3 - 6 * beta[j] * pfqn_ca(gammatilde,m)
            m = zeros((1,p))
            m[j] = 5
            A3 = A3 - 20 * beta[j] ** 2 * pfqn_ca(gammatilde,m)
            m = zeros((1,p))
            m[j] = 6
            A3 = A3 - 15 * beta[j] ** 3 * pfqn_ca(gammatilde,m)
            for k in setdiff1d(arange(0,p+1),j):
                m = zeros((1,p))
                m[j] = 4
                m[k] = 2
                A3 = A3 - 2 * beta[j] * beta[k] * pfqn_ca(gammatilde,m)
                m = zeros((1,p))
                m[j] = 2
                m[k] = 3
                A3 = A3 - 3 * beta[j] ** 2 * beta[k] * pfqn_ca(gammatilde,m)
                for l in setdiff1d(arange(0,p+1),array([j,k])):
                    m = zeros((1,p))
                    m[j] = 2
                    m[k] = 2
                    m[l] = 2
                    A3 = A3 - (1 / 6) * beta[j] * beta[k] * beta(l) * pfqn_ca(gammatilde,m)

    I = array([A0,A1 / N,A2 / N ** 2])

    lGn = - sum(factln(K)) + sum(multiply(K,log(sum(Z, 0)))) + log(sum(I)) - sum(log(alpha))
    Gn = exp(lGn)
    if not isfinite(lGn) :
        Gn = NaN
        lGn = NaN

    return Gn,lGn