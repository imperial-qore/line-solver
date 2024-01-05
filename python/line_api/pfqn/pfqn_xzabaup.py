from numpy import *

def pfqn_xzabaup(L = None,N = None,Z = None):
    XN = amin(array([1 / amax(L),N / (sum(L) + Z)]))
    return XN

