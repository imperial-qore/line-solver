import numpy as np
from numpy import *
from ..util import *

def pfqn_nc_sanitize(lambda_=None, L=None, N=None, Z=None):
    Tol = 1e-14
    # erase empty classes
    nnzclasses = nonzero(N)
    L = L[:, nnzclasses]
    N = N[nnzclasses]
    Z = Z[nnzclasses]
    lambda_ = lambda_[:, nnzclasses]
    # erase ill-defined classes
    zeroclasses = nonzero(np.isclose(L[:, nnzclasses] + Z[:, nnzclasses]))
    L[:, zeroclasses] = []
    N[:, zeroclasses] = []
    Z[:, zeroclasses] = []
    lambda_[:, zeroclasses] = []

    lGremaind = 0
    # find zero demand classes
    zerodemands = find(L < Tol)
    if not len(zerodemands) == 0:
        lGremaind = lGremaind + N(zerodemands) * transpose(log(Z(zerodemands))) - sum(log(N(zerodemands)))
        L[:, zerodemands] = []
        Z[:, zerodemands] = []
        N[:, zerodemands] = []

    # rescale demands
    Lmax = amax(L, [], 1)

    L = L / matlib.repmat(Lmax, L.shape[0], 1)
    Z = Z / matlib.repmat(Lmax, Z.shape[0], 1)
    lGremaind = lGremaind + N * transpose(log(Lmax))
    # sort from smallest to largest think time
    L = L[Z.argsort()]
    N = N[Z.argsort()]
    Z = Z[Z.argsort()]
    # ensure zero think time classes are anyway first
    zerothinktimes = find(Z < Tol)
    nonzerothinktimes = setdiff1d(arange(0, L.shape[1] + 1), zerothinktimes)
    L = L[:, array([zerothinktimes, nonzerothinktimes])]
    N = N[:, array([zerothinktimes, nonzerothinktimes])]
    Z = Z[:, array([zerothinktimes, nonzerothinktimes])]
    return lambda_, L, N, Z, lGremaind
