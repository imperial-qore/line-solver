from numpy import *

def oner(N = None,r = None):
    # N=ONER(N,r)
    # Decrement element in position of r of input vector

    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    for s in array([r]):
        N[s] = N[s] - 1

    return N
