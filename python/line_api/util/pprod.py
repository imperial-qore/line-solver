import numpy as np
from numpy import *

def pprod(n = None,N = None):
    # [N]=PPROD(N,N)
    # sequentially generate all vectors n: 0<=n<=N
    # n=pprod(N) - init
    # n=pprod(n,N) - next state
    if N is None:
        N = n
        n = zeros(N.size, dtype=np.uint32)
        return n

    R = len(N)
    if sum(n == N) == R:
        n = - 1
        return n

    s = R
    while s >= 0 and n[s-1] == N[s-1]:
        n[s-1] = 0
        s = s - 1

    if s == 0:
        return n

    n[s-1] = n[s-1] + 1
    return n
