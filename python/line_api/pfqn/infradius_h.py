from numpy import *
from numpy.matlib import repmat
from . import pfqn_gld

def infradius_h(x = None,L = None,N = None,alpha = None):
    M = L.shape[0]
    Nt = sum(N)
    beta = N / Nt
    t = exp(x) / (1 + exp(x))
    tb = sum(multiply(beta,exp(x)) / (1 + exp(x)), 1)
    h = lambda x = None: pfqn_gld(sum(multiply(L,repmat(exp(2 * pi * 1j * (t - tb)),M,1)), 1),Nt,alpha) * prod(exp(x) / (1 + exp(x)) ** 2, 1)
    y = zeros((x.shape[0],1))
    for i in arange(0,x.shape[0]+1):
        y[i] = real(h(x[i,:]))

    return y