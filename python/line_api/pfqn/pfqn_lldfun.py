from numpy import *
def pfqn_lldfun(n = None,lldscaling = None,nservers = None):
    # R = PFQN_LLDFUN(N,MU,C)

    # AMVA-QD queue-dependence function

    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.

    M = len(n)
    r = ones((M,1))
    smax = lldscaling.shape[1]
    alpha = 20

    for i in arange(0,M+1):
        ## handle servers
        if nservers is not None:
            if isinf(nservers[i]):
                r[i] = 1
            else:
                r[i] = r[i] / softmin(n[i],nservers[i],alpha)
                if isnan(r[i]):
                    r[i] = 1 / amin(n[i],nservers[i])
        ## handle generic lld
        if not len(lldscaling)==0  and range_(lldscaling[i,:]) > 0:
        r[i] = r[i] / interp1(arange(0,smax+1),lldscaling(i,arange(0,smax+1)),n[i],'spline')

return r