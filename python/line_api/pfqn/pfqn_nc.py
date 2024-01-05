from numpy import *
from ..util import *

def pfqn_nc(lambda_ = None,L = None,N = None,Z = None,varargin = None):
    # [LG,X,Q] = PFQN_NC(L,N,Z,VARARGIN)

    # L: Service demand matrix
    # N: Population vector
    # Z: Think time vector
    # varargin: solver options (e.g., SolverNC.defaultOptions)

    # LG: Logarithm of normalizing constant
    # X: System throughputs
    # Q: Mean queue-lengths

    options = Solver.parseOptions(varargin,SolverNC.defaultOptions)
    # backup initial parameters
    Rin = len(N)
    if sum(N) == 0 or len(N)==0:
        lG = 0
        X = []
        Q = []
        return lG,X,Q

    if len(lambda_)==0:
        lambda_ = 0 * N

    X = []
    Q = []
    # compute open class contributions
    Qopen = []
    lGopen = 0
    for i in arange(0,L.shape[0]+1):
        Ut[i] = (1 - lambda_ * transpose(L[i,:]))
        if isnan(Ut[i]):
            Ut[i] = 0
        L[i,:] = L[i,:] / Ut[i]
        Qopen[i,:] = multiply(lambda_,L[i,:]) / Ut[i]
        #lGopen = lGopen + log(Ut[i]);

    Qopen[isnan[Qopen]] = 0
    # then erase open classes
    N[isinf[N]] = 0
    # first remove empty classes
    nnzClasses = find(N)
    lambda_ = lambda_([:,:,nnzClasses])
    L = L[:,nnzclasses]
    N = N[:,nnzclasses]
    Z = Z[:,nnzclasses]
    ocl = find(isinf(N))
    # then scale demands in [0,1], importat that stays before the other
    # simplications in case both D and Z are all very small or very large in a
    # given class, in which case the may look to filter but not if all of them
    # are at the same scale
    R = len(N)
    scalevec = ones((1,R))
    #switch options.method
    #    case {'adaptive','comom','default'}
    #        # no-op
    #    otherwise
    #end
    for r in arange(0,R+1):
        scalevec[r] = amax(array([[L[:,r]],[Z[:,r]]]))

        #end
        L = L / matlib.repmat(scalevec,L.shape[0],1)
        Z = Z / scalevec
        # remove stations with no demand
        Lsum = sum(L, 1)
        Lmax = amax(L,[],2)
        demStations = find((Lmax / Lsum) > Distrib.Zero)
        noDemStations = setdiff1d(arange(0,L.shape[0]+1),demStations)
        L = L(demStations,:)
        if any(N((sum(L, 0) + sum(Z, 0)) == 0) > 0):
            line_warning(mfilename,'The model has no positive demands in any class.')
            if len(Z)==0 or sum(Z) < options.tol:
                lG = 0
            else:
                lG = - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0)))) + N * transpose(log(scalevec))
            return lG,X,Q

        # update M and R
            if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
        # return immediately if degenerate case
        if len(L)==0 or sum(L) < options.tol:
            if len(Z)==0 or sum(Z) < options.tol:
                lG = lGopen
            else:
                lG = lGopen - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0)))) + N * transpose(log(scalevec))
            return lG,X,Q
        else:
            if M == 1 and (len(Z)==0 or sum(Z) < options.tol):
                lG = factln(sum(N)) - sum(factln(N)) + sum(multiply(N,log(sum(L, 0)))) + N * transpose(log(scalevec))
                return lG,X,Q

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
    Zz = Z(:,zeroDemandClasses)
    Z = Z(:,nonzeroDemandClasses)
    scalevecz = scalevec(nonzeroDemandClasses)
    # compute G for classes No with non-zero demand
    lGnzdem,Xnnzdem,Qnnzdem = compute_norm_const(L,N,Z,options)
    if len(Xnnzdem)==0:
        X = []
        Q = []
    else:
        zClasses = setdiff1d(arange(0,Rin+1),nnzClasses)
        Xz = zeros((1,len(zClasses)))
        Xnnz = zeros((1,len(nnzClasses)))
        Xnnz[zeroDemandClasses] = Nz / sum(Zz, 0) / scalevec(zeroDemandClasses)
        Xnnz[nonzeroDemandClasses] = Xnnzdem / scalevec(nonzeroDemandClasses)
        X[1,array[[zClasses,nnzClasses]]] = array([Xz,Xnnz])
        X[ocl] = lambda_(ocl)
        Qz = zeros((Qnnzdem.shape[0],len(zClasses)))
        Qnnz = zeros((Qnnzdem.shape[0],len(nnzClasses)))
        Qnnz[:,zeroDemandClasses] = 0
        Qnnz[:,nonzeroDemandClasses] = Qnnzdem
        Q[noDemStations,:] = 0
        Q[demStations,array[[zClasses,nnzClasses]]] = array([Qz,Qnnz])
        Q[:,ocl] = Qopen(:,ocl)

        # scale back to original demands
        lG = lGopen + lGnzdem + lGzdem + N * transpose(log(scalevecz))
        return lG,X,Q


    def compute_norm_const(L = None,N = None,Z = None,options = None):
        # LG = COMPUTE_NORM_CONST(L,N,Z,OPTIONS)
        # Auxiliary script that computes LG after the initial filtering of L,N,Z

            if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
        X = []
        Q = []
        if array(['ca']) == options.method:
            __,lG = pfqn_ca(L,N,sum(Z, 0))
        else:
            if array(['adaptive','default']) == options.method:
                if M > 1:
                    if R == 1 or (R <= 3 and sum(N) < 50):
                        __,__,__,__,lG = pfqn_mva(L,N,sum(Z, 0))
                    else:
                        if M > R:
                            __,lG = pfqn_kt(L,N,sum(Z, 0))
                        else:
                            __,lG = pfqn_le(L,N,sum(Z, 0))
                else:
                    if sum(Z, 0) == 0:
                        lG = - N * transpose(log(L))
                    else:
                        if N < 10000:
                            #[lG] = pfqn_comomrm(L,N,sum(Z,1));
                            __,lG = pfqn_mmint2_gausslegendre(L,N,sum(Z, 0))
                        else:
                            __,lG = pfqn_le(L,N,sum(Z, 0))
            else:
                if array(['sampling']) == options.method:
                    if M == 1:
                        __,lG = pfqn_mmsample2(L,N,sum(Z, 0),options.samples)
                    else:
                        if M > R:
                            __,lG = pfqn_mci(L,N,sum(Z, 0),options.samples,'imci')
                        else:
                            __,lG = pfqn_ls(L,N,sum(Z, 0),options.samples)
                else:
                    if array(['mmint2']) == options.method:
                        if L.shape[0] > 1:
                            line_error(mfilename,sprintf('The %s method requires a model with a delay and a single queueing station.',options.method))
                        __,lG = pfqn_mmint2_gausslegendre(L,N,sum(Z, 0))
                    else:
                        if 'cub' == options.method:
                            __,lG = pfqn_cub(L,N,sum(Z, 0))
                        else:
                            if 'kt' == options.method:
                                __,lG = pfqn_kt(L,N,sum(Z, 0))
                            else:
                                if 'le' == options.method:
                                    __,lG = pfqn_le(L,N,sum(Z, 0))
                                else:
                                    if 'ls' == options.method:
                                        __,lG = pfqn_ls(L,N,sum(Z, 0),options.samples)
                                    else:
                                        if array(['mci','imci']) == options.method:
                                            __,lG = pfqn_mci(L,N,sum(Z, 0),options.samples,'imci')
                                        else:
                                            if array(['mva']) == options.method:
                                                __,__,__,__,lG = pfqn_mva(L,N,sum(Z, 0))
                                            else:
                                                if 'mom' == options.method:
                                                    if len(N) > 1:
                                                        try:
                                                            __,lG,X,Q = pfqn_mom(L,N,Z)
                                                        finally:
                                                            pass
                                                    else:
                                                        X,Q,__,__,lG = pfqn_mva(L,N,Z)
                                                    #case {'exact'}
                                                    #if M>=R || sum(N)>20 || sum(Z)>0
                                                    #    [~,lG] = pfqn_ca(L,N,sum(Z,1));
                                                    #else
                                                    #    [~,lG] = pfqn_recal(L,N,sum(Z,1));# implemented with Z=0
                                                    #end
                                                else:
                                                    if array(['exact','comom']) == options.method:
                                                        if R > 1:
                                                            try:
                                                                # comom has a bug in computing X, sometimes the
                                                                # order is switched
                                                                if M > 1:
                                                                    __,lG,X,Q = pfqn_comombtf_java(L,N,Z)
                                                                else:
                                                                    lG = pfqn_comomrm(L,N,Z)
                                                            finally:
                                                                pass
                                                        else:
                                                            X,Q,__,__,lG = pfqn_mva(L,N,Z)
                                                    else:
                                                        if array(['pana','panacea','pnc']) == options.method:
                                                            __,lG = pfqn_panacea(L,N,sum(Z, 0))
                                                            if isnan(lG):
                                                                line_warning(mfilename,'Model is not in normal usage, panacea cannot continue.')
                                                        else:
                                                            if 'propfair' == options.method:
                                                                __,lG = pfqn_propfair(L,N,sum(Z, 0))
                                                            else:
                                                                if array(['recal']) == options.method:
                                                                    if sum(Z) > 0:
                                                                        line_error(mfilename,'RECAL is currently available only for models with non-zero think times.')
                                                                    __,lG = pfqn_recal(L,N,sum(Z, 0))
                                                                else:
                                                                    if 'rgf' == options.method:
                                                                        if sum(Z) > 0:
                                                                            line_error(mfilename,'RGF is defined only for models with non-zero think times.')
                                                                        __,lG = pfqn_rgf(L,N)
                                                                    else:
                                                                        line_error(mfilename,sprintf('Unrecognized method: %s',options.method))

        return lG,X,Q