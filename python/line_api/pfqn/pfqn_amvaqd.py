from numpy import *
from ..util import *
def pfqn_amvaqd(L = None,N = None,ga = None,be = None,Q0 = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    Q = zeros((M,R))
    if len(varargin) < 5:
        Q = multiply(L / matlib.repmat(sum(L, 0),M,1),matlib.repmat(N,M,1))
    else:
        Q = Q0

    delta = (sum(N) - 1) / sum(N)
    deltar = (N - 1) / N
    Q_1 = Q * 10
    tol = 1e-06
    iter = 0
    while amax(amax(abs(Q - Q_1))) > tol:

        iter = iter + 1
        Q_1 = Q
        for k in arange(0,M+1):
            for r in arange(0,R+1):
                Ak[r][k,1] = 1 + delta * sum(Q[k,:])
                Akr[k,r] = 1 + deltar[r] * Q[k,r]
        #    Q
        for r in arange(0,R+1):
            g = ga(Ak[r])
            b = be(Akr)
            for k in arange(0,M+1):
                C[k,r] = L[k,r] * g[k] * b[k,r] * (1 + delta * sum(Q[k,:]))
                X[r] = N[r] / sum(C[:,r])
                for k in arange(0,M+1):
                    Q[k,r] = X[r] * C[k,r]
                    U[k,r] = L[k,r] * g[k] * b[k,r] * X[r]


        return Q,X,U,iter

    def pfqn_qdamva_iter(sn = None,Qchain = None,Xchain = None,Uchain = None,STchain = None,Vchain = None,Nchain = None,SCVchain = None,options = None):
        M = sn.nstations
        K = sn.nchains
        nservers = sn.nservers
        schedparam = sn.schedparam
        lldscaling = sn.lldscaling
        cdscaling = sn.cdscaling
        Wchain = zeros((M,K))
        Uhiprio = zeros((M,K))

        Nt = sum(Nchain(isfinite(Nchain)))
        delta = (Nt - 1) / Nt
        deltaclass = (Nchain - 1) / Nchain
        deltaclass[isinf[Nchain]] = 1
        nnzclasses = find(Nchain > 0)
        ocl = find(isinf(Nchain))
        ccl = find(logical_and(isfinite(Nchain),Nchain) > 0)
        ## evaluate lld and cd correction factors
        totArvlQlenSeenByOpen = zeros((M,1))
        interpTotArvlQlen = zeros((M,1))
        totArvlQlenSeenByClosed = zeros((M,K))
        selfArvlQlenSeenByOpen = zeros((M,K))
        selfArvlQlenSeenByClosed = zeros((M,K))
        for k in arange(0,M+1):
            totArvlQlenSeenByOpen[k] = sum(Qchain(k,nnzclasses))
            interpTotArvlQlen[k] = delta * sum(Qchain(k,nnzclasses))
            for r in nnzclasses:
                totArvlQlenSeenByClosed[k,r] = deltaclass[r] * Qchain[k,r] + sum(Qchain(k,setdiff1d(nnzclasses,r)))
                selfArvlQlenSeenByOpen[k,r] = Qchain[k,r]
                selfArvlQlenSeenByClosed[k,r] = deltaclass[r] * Qchain[k,r]

        if not len(lldscaling)==0  or not len(cdscaling)==0 :
            qdterm = pfqn_qdfun(1 + interpTotArvlQlen,lldscaling,nservers)
            for r in nnzclasses:
                cdterm = ones((M,K))
                if isfinite(NchaiN[r]):
                    cdterm[:,r] = pfqn_cdfun(1 + selfArvlQlenSeenByClosed,cdscaling)
                else:
                    cdterm[:,r] = pfqn_cdfun(1 + selfArvlQlenSeenByOpen,cdscaling)

        STeff = STchain

        Wchain = zeros((M,K))
        ## compute response times from current queue-lengths
        for r in nnzclasses:
            sd = setdiff1d(nnzclasses,r)
            for k in arange(0,M+1):
                if not len(lldscaling)==0  or not len(cdscaling)==0 :
                    STeff[k,r] = STchain[k,r] * qdterm[k] * cdterm[k,r]
                if SchedStrategy.ID_INF == sn.schedid[k]:
                    Wchain[k,r] = STeff[k,r]
                else:
                    if SchedStrategy.ID_PS == sn.schedid[k]:
                        if array(['default','amva','amva.qd','qd']) == options.method:
                            Wchain[k,r] = Wchain[k,r] + STeff[k,r] * (1 + interpTotArvlQlen[k])
                        else:
                            if array(['amva.bs','bs']) == options.method:
                                Wchain[k,r] = STeff[k,r] * (nservers[k] - 1) / nservers[k]
                                if ismember(r,ocl):
                                    Wchain[k,r] = Wchain[k,r] + STeff[k,r] * (1 / nservers[k]) * (1 + totArvlQlenSeenByOpen[k])
                                else:
                                    Wchain[k,r] = Wchain[k,r] + STeff[k,r] * (1 / nservers[k]) * (1 + totArvlQlenSeenByClosed[k])
                    else:
                        if SchedStrategy.ID_DPS == sn.schedid[k]:
                            weight = schedparam[k,:]
                            if nservers[k] > 1:
                                line_error(mfilename,'Multi-server DPS not supported yet in AMVA solver.')
                            else:
                                tss = Inf
                                Uhiprio[k,r] = sum(Uchain(k,weight[k,:] > tss * weight[k,r]))
                                STeff[k,r] = STeff[k,r] / (1 - Uhiprio[k,r])
                                Wchain[k,r] = STeff[k,r] * selfArvlQlenSeenByClosed[k,r]
                                for s in setdiff1d(sd,r):
                                    if weight[s] == weight[r]:
                                        Wchain[k,r] = Wchain[k,r] + STeff[k,r] * Qchain(k,s)
                                    else:
                                        if weight[s] / weight[r] < tss:
                                            Wchain[k,r] = Wchain[k,r] + STeff[k,r] * Qchain(k,s) * weight[s] / weight[r]
                                        else:
                                            if weight[s] / weight[r] > tss:
                                                # if there is time-scale separation, do nothing
                                                # all is accounted for by 1/(1-Uhiprio)
                                                pass
                        else:
                            if array([SchedStrategy.ID_FCFS,SchedStrategy.ID_SIRO,SchedStrategy.ID_LCFSPR]) == sn.schedid[k]:
                                if STeff[k,r] > 0:
                                    if nservers[k] > 1:
                                        B = ((multiply(multiply(multiply(deltaclass,Xchain),Vchain[k,:]),STeff[k,:])) / nservers[k])
                                        else:
                                        B = ones((1,K))
                                    if nservers[k] == 1 and (not len(lldscaling)==0  or not len(cdscaling)==0 ):
                                        Wchain[k,r] = STeff[k,r] * (1 + SCVchain[k,r]) / 2
                                        Wchain[k,r] = Wchain[k,r] + (STeff[k,r] * deltaclass[r] * Qchain[k,r] + STeff(k,sd) * transpose(Qchain(k,sd)))
                                    else:
                                        Wchain[k,r] = STeff[k,r] * (nservers[k] - 1) / nservers[k]
                                        Wchain[k,r] = Wchain[k,r] + (1 / nservers[k]) * STeff[k,r] * (1 + SCVchain[k,r]) / 2
                                        Wchain[k,r] = Wchain[k,r] + (1 / nservers[k]) * (STeff[k,r] * deltaclass[r] * Qchain[k,r] / B[r] + multiply((STeff(k,sd) / B(sd)),transpose(Qchain(k,sd))))

        return Wchain