
from numpy import *
from scipy.special import factorial

from ..util import factln


def pfqn_grnmol(L = None,N = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    G = 0
    S = ceil(sum(N) - 1) / 2
    H = zeros((1 + S,1))
    for i in arange(0,S+1):
        c[1 + i] = 2 * (S - i) + M
        w[1 + i] = 2 ** - (2 * S) * (- 1) ** i * c(1 + i) ** (2 * S + 1) / factorial[i] / factorial(i + c(1 + i))
        s,bvec,SD,D = sprod(M,S - i)
        bvec = transpose(bvec)
        while bvec[1] >= 0:

            H[1 + i] = H[1 + i] + prod((((2 * bvec + 1) / c(1 + i)) * L) ** N)
            s,bvec = sprod(s,SD,D)
            bvec = transpose(bvec)

        G = G + w(1 + i) * H[1 + i]

    G = G * exp(factln(sum(N) + M - 1) - sum(factln(N)))
    return G