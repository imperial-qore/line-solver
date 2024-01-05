from numpy import *
from ..util import *

def pfqn_mmint2_gausslaguerre(L = None,N = None,Z = None,m = None):
    # [G,LOGG] = PFQN_MMINT2_GAUSSLAGUERRE(L,N,Z,m)

    # Integrate with Gauss-Laguerre

    if len(varargin) < 4:
        m = 1



    if len(gausslaguerreNodes)==0:
        gausslaguerreNodes,gausslaguerreWeights = gengausslegquadrule(300,10 ** - 5)

    lambda_ = 0 * N
    nonzeroClasses = find(N)
    # repairmen integration
    f = lambda u = None: N(nonzeroClasses) * transpose(log(Z(nonzeroClasses) + L(nonzeroClasses) * u))
    x = gausslaguerreNodes
    w = gausslaguerreWeights
    n = amin(300,2 * sum(N) + 1)
    F = zeros((x.shape,x.shape))
    for i in arange(0,len(x)+1):
        F[i] = (m - 1) * log(x[i]) + f(x[i])

    g = log(w) + F - sum(factln(N)) - factln(m - 1)
    lG = log(sum(exp(g)))
    if not isfinite(lG) :
        lG = logsumexp(g)

    G = exp(lG)
    return G,lG