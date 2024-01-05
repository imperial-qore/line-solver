from numpy import *

def pfqn_mvams(lambda_ = None,L = None,N = None,Z = None,mi = None,S = None):
    # [XN,QN,UN,CN,LOGG]=PFQN_MVAMS(LAMBDA,L,N,Z,MI,S)

    # this is a general purpose script to handle mixed qns with multi-server nodes
    # S[i] number of servers in station i
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]

    Ntot = sum(N(isfinite(N)))
    mu = ones((M,Ntot))
    if len(varargin) < 6:
        S = ones((M,1))

    if len(varargin) < 5:
        mi = ones((M,1))

    if len(Z)==0:
        Z = zeros((1,R))

    for i in arange(0,M+1):
        mu[i,:] = amin(arange(0,Ntot+1),S[i] * ones((1,Ntot)))

    if amax(S(isfinite(S))) == 1:
        if any(isinf(N)):
            XN,QN,UN,CN,lG = pfqn_mvamx(lambda_,L,N,Z,mi)
        else:
            XN,QN,UN,CN,lG = pfqn_mva(L,N,Z,mi)
    else:
        if any(isinf(N)):
            if amax(mi) == 1:
                lG = NaN
                XN,QN,UN,CN = pfqn_mvaldms(lambda_,L,N,Z,S)
            else:
                line_error(mfilename,'Queue replicas not available in exact MVA for mixed models.')
        else:
            XN,QN,UN,CN,lG = pfqn_mvald(L,N,Z,mu)
            lG = lG(end())

    return XN,QN,UN,CN,lG