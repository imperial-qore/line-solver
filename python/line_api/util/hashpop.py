from numpy import *

def hashpop(n = None, N = None,R = None, prods = None):
    # IDX=HASHPOP(N,N,R,PRODS)

    # hash a population vector in n: 0<=n<=N
    idx = 1
    if R is None:
        R = n.size
        for r in arange(0,R):
            idx = idx + prod(N[arange(0,r)] + 1) * n[r]
        return idx-1
    else:
        for r in arange(0,R):
            idx = idx + prods[r] * n[r]
    return idx-1
