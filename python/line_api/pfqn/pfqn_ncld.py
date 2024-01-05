from numpy import *
from ..util import *
def pfqn_ncld(L = None,N = None,Z = None,mu = None,varargin = None):
    # [LGN] = PFQN_NCLD(L,N,Z,VARARGIN)

    options = Solver.parseOptions(varargin,SolverNC.defaultOptions)
    # backup initial parameters

    X = []
    Q = []
    mu = mu(:,arange(0,sum(N)+1))
    # first remove empty classes
    nnzClasses = find(N)
    L = L[:,nnzclasses]
    N = N[:,nnzclasses]
    Z = Z[:,nnzclasses]
    # then scale demands in [0,1], importat that stays before the other
    # simplications in case both D and Z are all very small or very large in a
    # given class, in which case the may look to filter but not if all of them
    # are at the same scale
    R = len(N)
    scalevec = ones((1,R))
    for r in arange(0,R+1):
        scalevec[r] = amax(array([[L[:,r]],[Z[:,r]]]))

        L = L / matlib.repmat(scalevec,L.shape[0],1)
        Z = Z / scalevec
        # remove stations with no demand
        Lsum = sum(L, 1)
        Lmax = amax(L,[],2)
        demStations = find((Lmax / Lsum) > Distrib.Zero)
        L = L(demStations,:)
        mu = mu(demStations,:)
        # if there is a class with jobs but with L and Z all zero
        if any(N((sum(L, 0) + sum(Z, 0)) == 0) > 0):
            line_warning(mfilename,'The model has no positive demands in any class.')
            if len(Z)==0 or sum(Z) < options.tol:
                lG = 0
            else:
                lG = - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0)))) + N * transpose(log(scalevec))
            return lG,G

        # update M and R
            if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
        # return immediately if the model is a degenerate case
        if len(L)==0 or sum(L) < options.tol:
            if len(Z)==0 or sum(Z) < options.tol:
                lG = 0
            else:
                lG = - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0)))) + N * transpose(log(scalevec))
            return lG,G
        else:
            if M == 1 and (len(Z)==0 or sum(Z) < options.tol):
                lG = factln(sum(N)) - sum(factln(N)) + sum(multiply(N,log(sum(L, 0)))) + N * transpose(log(scalevec)) - sum(log(mu(arange(0,sum(N)+1))))
                return lG,G

        # determine contribution from jobs that permanently loop at delay
        zeroDemandClasses = find(sum(L, 0) < options.tol)

        nonzeroDemandClasses = setdiff1d(arange(0,R+1),zeroDemandClasses)
        if len(sum(Z, 0))==0 or all(sum(Z(:,zeroDemandClasses), 0) < options.tol):
        lGzdem = 0
        Nz = 0
    else:
        if len(zeroDemandClasses)==0:
            lGzdem = 0
            Nz = 0
        else:
            Nz = N(zeroDemandClasses)
            lGzdem = - sum(factln(Nz)) + sum(multiply(Nz,log(sum(Z(:,zeroDemandClasses), 0)))) + Nz * transpose(log(scalevec(zeroDemandClasses)))

    L = L(:,nonzeroDemandClasses)
    N = N(nonzeroDemandClasses)
    Z = Z(:,nonzeroDemandClasses)
    scalevecz = scalevec(nonzeroDemandClasses)
    # compute G for classes No with non-zero demand
    lGnnzdem = compute_norm_const_ld(L,N,Z,mu,options)
    # scale back to original demands
    lG = lGnnzdem + lGzdem + N * transpose(log(scalevecz))
    G = exp(lG)
    return lG,G

def compute_norm_const_ld(L = None,N = None,Z = None,mu = None,options = None):
    # LG = COMPUTE_NORM_CONST_LD(L,N,Z,OPTIONS)
    __,R = L.shape
    D = Z.shape[0]

    Lz = array([[L],[Z]])
    muz = array([[mu],[matlib.repmat(arange(0,mu.shape[1]+1),D,1)]])
    if array(['default','exact']) == options.method:
        if R == 1:
            lG = pfqn_gldsingle(Lz,N,muz,options)
        else:
            __,lG = pfqn_gld(Lz,N,muz,options)
    else:
        if 'rd' == options.method:
            lG = pfqn_rd(L,N,Z,mu,options)
        else:
            if array(['nrp','nr.probit']) == options.method:
                lG = pfqn_nrp(L,N,Z,mu,options)
            else:
                if array(['nrl','nr.logit']) == options.method:
                    lG = pfqn_nrl(L,N,Z,mu,options)
                else:
                    line_error(mfilename,sprintf('Unrecognized method for solving load-dependent models: %s',options.method))

    return lG