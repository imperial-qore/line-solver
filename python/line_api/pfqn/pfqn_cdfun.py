from numpy import *

def pfqn_cdfun(nvec=None, cdscaling=None):
    # R = PFQN_CDFUN(NVEC,PHI)
    #
    # AMVA-QD class-dependence function
    #
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    M = nvec.shape[0]
    r = ones((M, 1))
    if cdscaling is not None and not len(cdscaling) == 0:
        for i in arange(0, M):
            r[i] = 1 / cdscaling[i](nvec[i, :])

        return r
