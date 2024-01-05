from numpy import *
from ..util import *
def pfqn_mmsample2(L = None,N = None,Z = None,samples = None):
    # [G,LOGG,LOGGR] = PFQN_MMSAMPLE2(L,N,Z,SAMPLES)

    # Monte carlo sampling for normalizing constant of a repairmen model
    # based on McKenna-Mitra integral form
    R = len(N)
    # Scale so that all coefficients are >=1.
    scaleFactor = 1e-07 + amin(array([[L],[Z]]))
    L = L / scaleFactor
    Z = Z / scaleFactor
    nnzeros = arange(0,R+1)
    c = 0.5
    v = array([random.rand(1,ceil(c * samples)),logspace(0,5,ceil(samples * (1 - c)))])

    du = transpose(array([0,diff(v)]))
    u = matlib.repmat(transpose(v),1,len(nnzeros))
    ZL = log(multiply(matlib.repmat(Z + L(1,nnzeros),u.shape[0],1),u))
    lG = du + - transpose(v) + ZL * transpose(N)
    den = sum(factln(N))
    lG = amax(lG) - den

    lG = lG + sum(N) * log(scaleFactor)

    for r in arange(0,R+1):
        lGr = lG - ZL[:,r]
        lGr[r] = amax(lGr) - den + log(N[r]) + (sum(N) - 1) * log(scaleFactor)

    G = exp(lG)
    return G,lG,lGr
