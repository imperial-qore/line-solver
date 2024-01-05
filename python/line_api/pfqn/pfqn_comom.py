from numpy import *
from ..util import *

def pfqn_comom(L = None,N = None,Z = None):
    if len(varargin) < 3:
        Z = 0 * N

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    # rescale demands
    Lmax = L

    Lmax[Lmax < Distrib.Tol] = Z(Lmax < Distrib.Tol)

    Lmax = amax(Lmax,[],1)
    L = L / matlib.repmat(Lmax,M,1)
    Z = Z / matlib.repmat(Lmax,M,1)
    # sort from smallest to largest
    #[~,rsort] = sort(Z,'ascend');
    #L=L(:,rsort);
    #Z=Z(:,rsort);
    # prepare comom data structures
    Dn = multichoose(R,M)
    Dn[:,R] = 0
    Dn = sortbynnzpos(Dn)
    # initialize
    nvec = zeros((1,R))
    h = ginit(L)
    lh = log(h) + factln(sum(nvec) + M - 1) - sum(factln(nvec))
    h = exp(lh)
    scale = zeros((1,sum(N)))
    # iterate
    for r in arange(0,R+1):
        for Nr in arange(0,N[r]+1):
            nvec[r] = nvec[r] + 1
            if Nr == 1:
                A,B,DA = genmatrix(L,nvec,Z,r)
            else:
                A = A + DA
            b = B * h * nvec[r] / (sum(nvec) + M - 1)
            h = linalg.solve(A,b)
            nt = sum(nvec)
            scale[nt] = abs(sum(__builtint__.sorted(h)))
            h = abs(h) / scale(nt)

    # unscale and return the log of the normalizing constant
    lG = log(h(end() - (R - 1))) + factln(sum(N) + M - 1) - sum(factln(N)) + N * transpose(log(Lmax)) + sum(log(scale))
    return lG

def genmatrix(L = None,N = None,Z = None,r = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    A = zeros((nchoosek(M + R - 1,M) * (M + 1),nchoosek(M + R - 1,M) * (M + 1)))
    DA = zeros((nchoosek(M + R - 1,M) * (M + 1),nchoosek(M + R - 1,M) * (M + 1)))
    B = zeros((nchoosek(M + R - 1,M) * (M + 1),nchoosek(M + R - 1,M) * (M + 1)))
    row = 0
    rr = []

    lastnnz = 0
    for d in arange(0,len(Dn)+1):
        hnnz = hashnnz(Dn(d,:),R)
        if hnnz != lastnnz:
            lastnnz = hnnz

    for d in arange(0,len(Dn)+1):
        if sum(Dn(d,arange([r],R - 1+1))) > 0:
            # dummy rows for unused norm consts
            for k in arange(0,M+1):
                row = row + 1
                col = hash(N,N - Dn(d,:),k + 1)
                A[row,col] = 1
                if sum(Dn(d,arange((r + 1),R - 1+1))) > 0:
                    col = hash(N,N - Dn(d,:),k + 1)
                    B[row,col] = 1
                else:
                    er = zeros((1,R))
                    er[r] = 1
                    col = hash(N,N - Dn(d,:) + er,k + 1)
                    B[row,col] = 1
        else:
            if sum(Dn(d,arange(0,r+1))) < M:
                for k in arange(0,M+1):
                    # add CE
                    row = row + 1
                    A[row,hash[N,N - Dn[d,:],k + 1]] = 1
                    A[row,hash[N,N - Dn[d,:],0 + 1]] = - 1
                    for s in arange(0,r - 1+1):
                        A[row,hash[N,oner[N - Dn[d,:],s],k + 1]] = - L(k,s)
                    B[row,hash[N,N - Dn[d,:],k + 1]] = L[k,r]
                for s in arange(0,(r - 1)+1):
                    # add PC to A
                    row = row + 1
                    n = N - Dn(d,:)
                    A[row,hash[N,n,0 + 1]] = N[s]
                    A[row,hash[N,oner[n,s],0 + 1]] = - Z[s]
                    for k in arange(0,M+1):
                        A[row,hash[N,oner[n,s],k + 1]] = - L(k,s)
                    B[row,:] = 0

    #add PC of class R
    for d in arange(0,len(Dn)+1):
        if sum(Dn(d,arange([r],R - 1+1))) <= 0:
            row = row + 1
            n = N - Dn(d,:)
            A[row,hash[N,n,0 + 1]] = N[r]
            DA[row,hash[N,n,0 + 1]] = 1
            B[row,hash[N,n,0 + 1]] = Z[r]
            for k in arange(0,M+1):
                B[row,hash[N,n,k + 1]] = L[k,r]

    rr = unique(rr)
    return A,B,DA,rr


def hashnnz(dn = None,R = None):
    val = 0
    for t in arange(0,R+1):
        if dn(t) == 0:
            val = val + 2 ** (t - 1)

    return val


def hash(N = None,n = None,i = None):
    if i == 1:
        col = Dn.shape[0] * M + matchrow(Dn,N - n)
    else:
        col = (matchrow(Dn,N - n) - 1) * M + i - 1

    return col


def ginit(L = None):
    e1 = zeros((1,R))
    g = zeros((Dn.shape[0] * (M + 1),1))
    for i in arange(0,M+1):
        g[hash[N,N,i + 1]] = 1

    g = g
    return g