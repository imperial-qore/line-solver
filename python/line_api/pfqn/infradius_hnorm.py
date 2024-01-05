from numpy import *
from numpy.matlib import repmat
from scipy.stats import norm
from . import pfqn_gld

def infradius_hnorm(x = None,L = None,N = None,alpha = None):
    M = L.shape[0]
    MU = 0
    SIGMA = 1
    Nt = sum(N)
    beta = N / Nt
    t = norm.cdf(x,MU,SIGMA)
    tb = sum(multiply(beta,t), 1)
    h = lambda x = None: pfqn_gld(sum(multiply(L,repmat(exp(2 * pi * 1j * (t - tb)),M,1)), 1),Nt,alpha) * prod(norm.pdf(x,MU,SIGMA), 1)
    y = zeros((x.shape[0],1))
    for i in arange(0,x.shape[0]+1):
        y[i] = real(h(x[i,:]))

    return y
