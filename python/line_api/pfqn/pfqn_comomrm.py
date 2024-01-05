from numpy import *
from ..util import *

def pfqn_comomrm(L = None,N = None,Z = None,m = None):
    # comom for a finite repairment model
    if L.shape[0] != 1:
        line_error(mfilename,'The solver accepts at most a single queueing station.')

    if len(varargin) < 4:
        m = 1

    lambda_ = 0 * N
    __,L,N,Z,lG0 = pfqn_nc_sanitize(lambda_,L,N,Z)
    zerothinktimes = nonzero(Z < 1e-06)
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    # initialize
    nvec = zeros((1,R))
    if any(zerothinktimes):
        nvec[zerothinktimes] = N(:,zerothinktimes)
        lh = []
        # these are trivial models with a single queueing station with demands all equal to one and think time 0
        lh[end() + 1,1] = (factln(sum(nvec) + m + 1 - 1) - sum(factln(nvec)))
        for s in find(zerothinktimes):
            nvec_s = oner(nvec,s)
            lh[end() + 1,1] = (factln(sum(nvec_s) + m + 1 - 1) - sum(factln(nvec_s)))
        lh[end() + 1,1] = (factln(sum(nvec) + m - 1) - sum(factln(nvec)))
        for s in find(zerothinktimes):
            nvec_s = oner(nvec,s)
            lh[end() + 1,1] = (factln(sum(nvec_s) + m - 1) - sum(factln(nvec_s)))
    else:
        lh = zeros((2,1))

    h = exp(lh)
    if len(zerothinktimes) == R:
        lGbasis = log(h)
        lG = lG0 + log(h(end() - R))
        return lG,lGbasis
    else:
        scale = ones((1,sum(N)))
        nt = sum(nvec)
        h_1 = h
        #iterate
        for r in arange((len(zerothinktimes) + 1),R+1):
            for Nr in arange(0,N[r]+1):
                nvec[r] = nvec[r] + 1
                if Nr == 1:
                    if r > len(zerothinktimes) + 1:
                        hr = zeros((2 * r,1))
                        hr[arange[1,[r - 1]+1]] = h(arange(0,(r - 1)+1))
                        hr[arange[[r + 1],[2 * r - 1]+1]] = h(arange(((r - 1) + 1),2 * (r - 1)+1))
                        h = hr
                        # update scalings
                        h[r] = h_1(1) / scale(nt)
                        h[end()] = h_1((r - 1) + 1) / scale(nt)
                    # CE for G+
                    A12[1,1] = - 1
                    # Class-1..(R-1) PCs for G
                    for s in arange(0,(r - 1)+1):
                        A12[1 + s,1] = N[s]
                        A12[1 + s,1 + s] = - Z[s]
                    # Class-R PCs
                    B2r = array([m * L(1,r) * eye[r],Z[r] * eye[r]])
                    # explicit formula for inv(C)
                    iC = - eye[r] / m
                    iC[1,:] = - 1 / m
                    iC[1] = 1
                    # explicit formula for F1r
                    F1r = zeros((2 * r,2 * r))
                    F1r[1,1] = 1
                    # F2r by the definition
                    F2r = array([[- iC * A12 * B2r],[B2r]])
                h_1 = h
                h = (F1r + F2r / nvec[r]) * h_1
                nt = sum(nvec)
                scale[nt] = abs(sum(__builtint__.sorted(h)))
                h = abs(h) / scale(nt)
        # unscale and return the log of the normalizing constant
        lG = lG0 + log(h(end() - (R - 1))) + sum(log(scale))
        lGbasis = log(h) + sum(log(scale))

    return lG,lGbasis
        