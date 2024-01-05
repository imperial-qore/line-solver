from numpy import *
from ..util import *
def pfqn_mmint2_gausslegendre(L = None,N = None,Z = None,m = None):
    # [G,LOGG] = PFQN_MMINT2_GAUSSLEGENDRE(L,N,Z,m)

    # Integrate McKenna-Mitra integral form with Gauss-Legendre in [0,1e6]
    if len(varargin) < 4:
        m = 1

    # nodes and weights generated with tridiagonal eigenvalues method in
    # high-precision using Julia:

    # using LinearAlgebra

    # function gauss(a, b, N)
    #    λ, Q = eigen(SymTridiagonal(zeros(N), [n / sqrt(4n^2 - 1) for n = 1:N-1]))
    #    @. (λ + 1) * (b - a) / 2 + a, [2Q[1, i]^2 for i = 1:N] * (b - a) / 2
    # end

    if len(gausslegendreNodes)==0:
        gausslegendreNodes = scipy.io.loadmat(which('gausslegendre-nodes.txt'))
        gausslegendreWeights = scipy.io.loadmat(which('gausslegendre-weights.txt'))

    # use at least 300 points
    n = amax(300,amin(len(gausslegendreNodes),2 * (sum(N) + m - 1) - 1))
    y = zeros((1,n))
    for i in arange(0,n+1):
        y[i] = N * transpose(log(Z + L * gausslegendreNodes[i]))

    g = log(gausslegendreWeights(arange(0,n+1))) - gausslegendreNodes(arange(0,n+1)) + y
    coeff = - sum(factln(N)) - factln(m - 1) + (m - 1) * sum(log(gausslegendreNodes(arange(0,n+1))))
    lG = log(sum(exp(g))) + coeff
    if not isfinite(lG) :
        lG = logsumexp(g) + coeff

    G = exp(lG)
    return G,lG