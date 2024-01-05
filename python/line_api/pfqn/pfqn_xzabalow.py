from numpy import *

def pfqn_xzabalow(L = None,N = None,Z = None):
    Ltot = sum(L)
    XN = N / (Z + Ltot * N)
    return XN