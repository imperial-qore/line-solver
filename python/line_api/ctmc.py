import numpy as np
    
def ctmc_makeinfgen(Q = None): 
    # Q=CTMC_MAKEINFGEN(Q)
    
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    A = Q - diag(diag(Q))
    Q = A - diag(np.sum(A, 2-1))
    return Q
    
    
def ctmc_rand(n = None): 
    # Q=CTMC_RAND(N)
    
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    Q = ctmc_makeinfgen(np.random.rand(n))
    return Q
    
def ctmc_simulate(Q = None,pi0 = None,n = None): 
    if len(pi0)==0:
        r = np.random.rand(len(Q),1)
        r = r / sum[r]
    
    __,st = np.amin(np.abs(rand - cumsum(pi0)))
    F = cumsum(Q - diag(diag(Q)),2)
    F = F / np.matlib.repmat(F(:,end()),1,len(F))
    for i in np.arange(1,n+1).reshape(-1):
        sts[i] = st
        soujt[i] = exprnd(- 1 / Q(st,st))
        st = 1 + np.amax(np.array([0,find(rand - F(st,:) > 0)]))
    
    return soujt,sts
    
def ctmc_ssg(sn = None,options = None): 
    stateSpace,stateSpaceHashed,qnc = State.spaceGenerator(sn,options.cutoff,options)
    nodeStateSpace = qnc.space
    sn.space = nodeStateSpace
    if options.verbose:
        line_printf('\nCTMC state space size: %d states. ',stateSpace.shape[1-1])
    
    if not isfield(options,'hide_immediate') :
        options.hide_immediate = True
    
    nstateful = sn.nstateful
    nclasses = sn.nclasses
    sync = sn.sync
    A = len(sync)
    stateSpaceAggr = np.zeros((stateSpaceHashed.shape,stateSpaceHashed.shape))
    # for all synchronizations
    for a in np.arange(1,A+1).reshape(-1):
        stateCell = cell(nstateful,1)
        for s in np.arange(1,stateSpaceHashed.shape[1-1]+1).reshape(-1):
            state = stateSpaceHashed(s,:)
            # update state cell array and SSq
            for ind in np.arange(1,sn.nnodes+1).reshape(-1):
                if sn.isstateful(ind):
                    isf = sn.nodeToStateful(ind)
                    stateCell[isf] = sn.space[isf](state(isf),:)
                    if sn.isstation(ind):
                        ist = sn.nodeToStation(ind)
                        __,nir = State.toMarginal(sn,ind,stateCell[isf])
                        stateSpaceAggr[s,np.arange[[[ist - 1] * nclasses + 1],ist * nclasses+1]] = nir
    
    return stateSpace,stateSpaceAggr,stateSpaceHashed,nodeStateSpace,sn
    
def ctmc_stochcomp(Q = None,I = None): 
    # [S,Q11,Q12,Q21,Q22,T] = CTMC_STOCHCOMP(Q,I)
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    if len(varargin) == 1:
        I = np.arange(1,np.ceil(len(Q) / 2)+1)
    
    Ic = setdiff(np.arange(1,len(Q)+1),I)
    Q11 = Q(I,I)
    Q12 = Q(I,Ic)
    Q21 = Q(Ic,I)
    Q22 = Q(Ic,Ic)
    #I = eye(size(Q22));
    T = np.linalg.solve((- Q22),Q21)
    T = Q12 * T
    S = Q11 + T
    return S,Q11,Q12,Q21,Q22,T
    
def ctmc_transient(Q = None,pi0 = None,t0 = None,t1 = None): 
    # [PI,T]=CTMC_TRANSIENT(Q,PI0,T0,T1)
    
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    if len(varargin) == 2:
        t1 = pi0
        t0 = 0
        pi0 = np.ones((1,len(Q)))
        pi0 = pi0 / sum(pi0)
    
    if len(varargin) == 3:
        t1 = t0
        t0 = 0
    
    t,np.pi = ode23(ctmc_transientode,np.array([t0,t1]),pi0)

return pi,t
    
def ctmc_transientode(t = None,pi = None): 
    # DPIDT=CTMC_TRANSIENTODE(T,PI)
    
    np.pi = np.transpose(np.pi)
    dpidt = np.pi * Q
    dpidt = dpidt
    return dpidt
    
def ctmc_uniformization(pi0 = None,Q = None,t = None,tol = None,maxiter = None): 
    # [PI,KMAX]=CTMC_UNIFORMIZATION(PI0,Q,T,TOL,MAXITER)
    
    # Copyright (c) 2012-2022, Imperial College London
    # All rights reserved.
    if len(varargin) < 4:
        tol = 1e-12
    
    if len(varargin) < 5:
        maxiter = 100
    
    q = 1.1 * np.amax(np.abs(diag(Q)))
    Qs = speye(Q.shape) + sparse(Q) / q
    k = 0
    s = 1
    r = 1
    iter = 0
    kmax = 1
    while iter < maxiter:

        iter = iter + 1
        k = k + 1
        r = r * (q * t) / k
        s = s + r
        if (1 - np.exp(- q * t) * s) <= tol:
            kmax = k
            break

    
    np.pi = pi0 * (np.exp(- q * t))
    P = pi0
    ri = np.exp(- q * t)
    for j in np.arange(1,kmax+1).reshape(-1):
        P = P * Qs
        ri = ri * (q * t / j)
        np.pi = np.pi + ri * P
    
    return np.pi,kmax