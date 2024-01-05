from numpy import *
from ..util import *

def pfqn_rgf(L = None,N = None,mi = None):
    if len(varargin) == 2:
        L,__,J = unique(L,'rows')
        mi = []
        for i in arange(0,L.shape[0]+1):
            mi[1,i] = sum(J == i)

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    Nck = zeros((sum(N) + sum(mi),sum(N) + sum(mi)))
    Dnk = cell(1,sum(N))
    Gsingle = zeros((1,2 * M + 1))
    iset = pprod(ones((1,R - 1)) * (M - 1))
    G = 0
    Ir = cell(R,1)
    Jr = zeros((M,1))
    ctr = 0
    gctr = 0
    while iset != - 1:

        g = rgfaux(L,N,mi,iset)
        G = G + g
        ctr = ctr + 1
        iset = pprod(iset,ones((1,R - 1)) * (M - 1))


    lG = log(G)

def rgfaux(L = None,N = None,mi = None,iset = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    X = L
    G = 0
    ## determine the loadings at the different steps
    for r in arange(R,2+- 1,- 1):
        for l in arange(0,M+1):
            if l != (iset(r - 1) + 1):
                for p in arange(0,(r - 1)+1):
                    X[l,p] = (X(iset(r - 1) + 1,r) * X(l,p) - X(l,r) * X(iset(r - 1) + 1,p)) / (X(iset(r - 1) + 1,r) - X(l,r))

    ## determine the multiplicities at the different steps
    # for each class we keep track of the couple (j,k), where k is the
    # index in the summation of the combinations that sum to j-1
    jstate = zeros((1,R))
    kstate = zeros((1,R))
    jmax = zeros((1,R))

    kmax = zeros((1,R))

    mistate = zeros((M,R))

    # initialize the elements of the state vectors in class R
    mistate[:,R] = mi
    r = R
    g = zeros((1,R))
    c = zeros((1,R))
    Jr = zeros((M,1))
    # initialize DFS recursion on jstate and kstate
    while r > 1:

        jmax[r] = mistate(iset(r - 1) + 1,r)
        jstate[r] = 1
        I[r] = mchoose(M - 1,jstate(R) - 1)
        kmax[r] = I[r].shape[0]
        kstate[r] = 1
        if iset(r - 1) + 1 > 1:
            Jr[arange[1,[iset[r - 1] + 1 - 1]+1],1] = I[r](kstate[r],arange(0,(iset(r - 1) + 1 - 1)+1))
        Jr[iset[r - 1] + 1,1] = N[r] + 1 - jstate[r]
        if iset(r - 1) + 1 < M:
            Jr[arange[[iset[r - 1] + 1 + 1],M+1],1] = I[r](kstate[r],arange(iset(r - 1) + 1,(M - 1)+1))
        mistate[:,r - 1] = mistate[:,r] + Jr
        # update g
        g[r] = (- 1) ** (jstate[r] + 1) * nck(mistate(iset(r - 1) + 1,r) + N[r] - jstate[r],N[r]) * X(iset(r - 1) + 1,r) ** N[r]
        for l in arange(0,M+1):
            if l != iset(r - 1) + 1:
                g[r] = g[r] * (X(iset(r - 1) + 1,r) / (X(iset(r - 1) + 1,r) - X(l,r))) ** mistate(l,r)
        c[r] = 1
        for l in arange(0,M+1):
            if l != iset(r - 1) + 1:
                c[r] = c[r] * nck(mistate(l,r) + Jr(l) - 1,Jr(l)) * (X(l,r) / (X(iset(r - 1) + 1,r) - X(l,r))) ** Jr(l)
        r = r - 1


    r = 2
    while r < R + 1:

        while jstate[r] <= jmax[r]:

            # update g
            g[r] = (- 1) ** (jstate[r] + 1) * nck(mistate(iset(r - 1) + 1,r) + N[r] - jstate[r],N[r]) * X(iset(r - 1) + 1,r) ** N[r]
            for l in arange(0,M+1):
                if l != iset(r - 1) + 1:
                    g[r] = g[r] * (X(iset(r - 1) + 1,r) / (X(iset(r - 1) + 1,r) - X(l,r))) ** mistate(l,r)
            while kstate[r] <= kmax[r]:

                # update multiplicities of class r
                if iset(r - 1) + 1 > 1:
                    Jr[arange[1,[iset[r - 1] + 1 - 1]+1],1] = I[r](kstate[r],arange(0,(iset(r - 1) + 1 - 1)+1))
                Jr[iset[r - 1] + 1,1] = N[r] + 1 - jstate[r]
                if iset(r - 1) + 1 < M:
                    Jr[arange[[iset[r - 1] + 1 + 1],M+1],1] = I[r](kstate[r],arange(iset(r - 1) + 1,(M - 1)+1))
                mistate[:,r - 1] = mistate[:,r] + Jr
                # update c
                c[r] = 1
                for l in arange(0,M+1):
                    if l != iset(r - 1) + 1:
                        c[r] = c[r] * nck(mistate(l,r) + Jr(l) - 1,Jr(l)) * (X(l,r) / (X(iset(r - 1) + 1,r) - X(l,r))) ** Jr(l)
                if r == 2:
                    G = G + prod(g(arange(2,R+1))) * prod(c(arange(2,R+1))) * gsingle(X(:,1),N(1),transpose(mistate(:,1)))
                    kstate[r] = kstate[r] + 1
                else:
                    break

            if r == 2:
                jstate[r] = jstate[r] + 1
                I[r] = mchoose(M - 1,jstate[r] - 1)
                kmax[r] = I[r].shape[0]
                kstate[r] = 1
            else:
                if r > 2:
                    r = r - 1
                    jmax[r] = mistate(iset(r - 1) + 1,r)
                    kstate[r] = 1
                    jstate[r] = 1
                    I[r] = mchoose(M - 1,jstate[r] - 1)
                    kmax[r] = I[r].shape[0]

        # find the largest class that needs updating
        r = r + 1
        if r > R:
            break
        while (kmax[r] <= kstate[r] and jmax[r] <= jstate[r]):

            r = r + 1
            if r > R:
                break

        # if completed go to the next iset element
        if r > R:
            break
        # update class r state
        if kmax[r] > kstate[r]:
            kstate[r] = kstate[r] + 1
        else:
            jstate[r] = jstate[r] + 1
            kstate[r] = 1
            I[r] = mchoose(M - 1,jstate[r] - 1)
            kmax[r] = I[r].shape[0]
        jmax[arange[1,[r - 1]+1]] = 0
        jstate[arange[1,[r - 1]+1]] = 0
        kmax[arange[1,[r - 1]+1]] = 0
        kstate[arange[1,[r - 1]+1]] = 0
        mistate[:,arange[1,[r - 1]+1]] = 0



def mchoose(n = None,k = None):
    I = Dnk[1,k + 1]
    if len(I)==0:
        I = multichoose(n,k)
        Dnk[1,k + 1] = I

    return I


def nck(c = None,k = None):
    b = Nck(c + 1,k + 1)
    if b == 0:
        b = nchoosek(c,k)
        Nck[c + 1,k + 1] = b

    return b


def gsingle(L = None,n = None,mult = None):
    L = round(L * 1000000.0) / 1000000.0
    row = matchrow(Gsingle(:,arange(0,2 * M+1)),array([transpose(L),mult]))
    if row == - 1:
        gctr = gctr + 1
        __,__,__,__,lG = pfqn_mva(L,n,0 * n,mult)
        G = exp(lG)
        Gsingle[end() + 1,:] = array([transpose(L),mult,G])
        return G

    G = Gsingle(row,2 * M + 1)
    return G
