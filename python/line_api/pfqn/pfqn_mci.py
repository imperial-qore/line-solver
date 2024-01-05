from numpy import *
from ..util import *
def pfqn_mci(D = None,N = None,Z = None,I = None,variant = None):
    # [G,LG,LZ] = PFQN_MCI(D,N,Z,I,VARIANT)

    # Normalizing constant estimation via Monte Carlo Integration

    # Syntax:
    # [lG,lZ] = gmvamcint(D,N,Z,I,Iest)
    # Input:
    # D - demands (queues x classes)
    # N - populations (1 x classes)
    # Z - think times (1 x classes)
    # I - samples

    # Output:
    # lG - estimate of logG
    # lZ - individual random samples

    # Note: if the script returns a floating point range exception,
    # double(log(mean(exp(sym(lZ))))) provides a better estimate of lG, but it
    # is very time consuming due to the symbolic operations.

    # Implementation: Giuliano Casale (g.casale@imperial.ac.uk), 16-Aug-2013

    if len(varargin) < 3:
        Z = 0 * N

    if len(varargin) < 4:
        I = 100000.0

    if len(varargin) < 5:
        variant = 'imci'

    M,R = D.shape
    if len(D)==0 or sum(D) < 0.0001:
        lGn = - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0))))
        G = exp(lGn)
        lZ = []
        return G,lG,lZ

    #tput = N./(Z+sum(D)+max(D).*(sum(N)-1)); # balanced job bounds
    #tput = N./Z; # balanced job bounds

    ## IMCI
    if strcmpi(variant,'imci'):
        tput = pfqn_amvabs(D,N,Z)
        util = D * transpose(tput)
        gamma = transpose(amax(0.01,1 - util))
    else:
        if strcmpi(variant,'mci'):
            tput = pfqn_amvabs(D,N,Z)
            util = D * transpose(tput)
            ## Original MCI
            for i in arange(0,len(util)+1):
                if util[i] > 0.9:
                    gamma[i] = 1 / sqrt(amax(N))
                else:
                    gamma[i] = 1 - util[i]
        else:
            if strcmpi(variant,'rm'):
                tput = N / (sum(D, 0) + Z + amax(D,1) * (sum(N) - 1))
                util = D * transpose(tput)
                ## Original MCI
                for i in arange(0,len(util)+1):
                    if util[i] > 0.9:
                        gamma[i] = 1 / sqrt(amax(N))
                    else:
                        gamma[i] = 1 - util[i]

    try:
        for r in arange(0,R+1):
            logfact[r] = sum(log(arange(0,N[r]+1)))
        # uniform sampling
        if len(VL)==0 or VL.shape[0] < I or VL.shape[1] < M:
            VL = log(random.rand(I,M))
        V = multiply(matlib.repmat(- 1.0 / gamma,I,1),VL(arange(0,I+1),arange(0,M+1)))
        ZI = matlib.repmat(Z,I,1)
        # importance sampling
        lZ = - (ones((1,M)) - gamma) * transpose(V) - sum(log(gamma)) - sum(logfact) + N * transpose(log(V * D + ZI))
        lG = logmeanexp(lZ)
        if isinf(lG):
            #    line_warning(mfilename,'Floating-point range exception, Monte Carlo integration will return an approximation.');
            lG = amax(lZ)
        G = exp(lG)
    finally:
        pass

    return G,lG,lZ