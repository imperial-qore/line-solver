from numpy import *
from ..util import *

def pfqn_mvaldmx(lambda_ = None,D = None,N = None,Z = None,mu = None,S = None):
    # [XN,QN,UN,CN] = PFQN_MVALDMX(LAMBDA,D,N,Z,MU,S)

    if len(varargin) < 5:
        mu = ones((D.shape[0],sum(N(isfinite(N)))))
        S = ones((D.shape[0],1))

    if mu.shape[1] < sum(N(isfinite(N))):
        line_error(mfilename,'MVALDMX requires to specify the load-dependent rates with one job more than the maximum closed population.')

    if any(N(find(lambda_)) > logical_and(0,isfinite(N(find(lambda_))))):
        line_error(mfilename,'Arrival rate cannot be specified on closed classes.')

    M,R = D.shape
    openClasses = find(isinf(N))
    closedClasses = setdiff1d(arange(0,len(N)+1),openClasses)
    XN = zeros((1,R))
    UN = zeros((M,R))
    CN = zeros((M,R))
    QN = zeros((M,R))
    lGN = 0
    mu[:,end() + 1] = mu(:,end())

    EC,E,Eprime = pfqn_mvaldmx_ec(lambda_,D,mu)
    C = len(closedClasses)

    Dc = D(:,closedClasses)
    Nc = N(closedClasses)
    Zc = Z(closedClasses)
    prods = zeros((1,C))

    for r in arange(0,C+1):
        prods[r] = prod(Nc(arange(0,r - 1+1)) + 1)

    # Start at nc=(0,...,0)
    nvec = pprod(Nc)
    # Initialize Pc
    Pc = zeros((M,1 + sum(Nc),prod(1 + Nc)))
    x = zeros((C,prod(1 + Nc)))
    w = zeros((M,C,prod(1 + Nc)))
    for i in arange(0,M+1):
        Pc[i,1 + 0,hashpop[nvec,Nc,C,prods]] = 1.0

    u = zeros((M,C))
    # Population recursion
    while nvec >= 0:

        hnvec = hashpop(nvec,Nc,C,prods)
        nc = sum(nvec)
        for i in arange(0,M+1):
            for c in arange(0,C+1):
                if nvec(c) > 0:
                    hnvec_c = hashpop(oner(nvec,c),Nc,C,prods)
                    # Compute mean residence times
                    for n in arange(0,nc+1):
                        w[i,c,hnvec] = w(i,c,hnvec) + Dc(i,c) * n * EC(i,n) * Pc(i,1 + (n - 1),hnvec_c)
        # Compute tput
        for c in arange(0,C+1):
            x[c,hnvec] = nvec(c) / (Zc(c) + sum(w(arange(0,M+1),c,hnvec)))
        for i in arange(0,M+1):
            for n in arange(0,nc+1):
                for c in arange(0,C+1):
                    if nvec(c) > 0:
                        hnvec_c = hashpop(oner(nvec,c),Nc,C,prods)
                        Pc[i,1 + n,hnvec] = Pc(i,1 + n,hnvec) + Dc(i,c) * EC(i,n) * x(c,hnvec) * Pc(i,1 + (n - 1),hnvec_c)
            Pc[i,1 + 0,hnvec] = amax(eps,1 - sum(Pc(i,1 + (arange(0,nc+1)),hnvec)))
        # now compute the normalizing constant
        last_nnz = find(nvec > 0,1,'last')
        if sum(nvec(arange(0,last_nnz - 1+1))) == sum(Nc(arange(0,last_nnz - 1+1))) and sum(nvec(arange((last_nnz + 1),C+1))) == 0:
            logX = log(XN(last_nnz))
            if not len(logX)==0 :
                lGN = lGN - logX
        nvec = pprod(nvec,Nc)


    # compute performance indexes at Nc for closed classes
    hnvec = hashpop(Nc,Nc,C,prods)
    for c in arange(0,C+1):
        hnvec_c = hashpop(oner(Nc,c),Nc,C,prods)
        for i in arange(0,M+1):
            u[i,c] = 0
            for n in arange(0,sum(Nc)+1):
                u[i,c] = u(i,c) + Dc(i,c) * x(c,hnvec) * Eprime(i,1 + n - 1) / E(i,1 + n - 1) * Pc(i,1 + n - 1,hnvec_c)

    # Throughput
    XN[closedClasses] = x(arange(0,C+1),hnvec)
    # Utilization
    UN[arange[1,M+1],closedClasses] = u(arange(0,M+1),arange(0,C+1))
    # Response time
    CN[arange[1,M+1],closedClasses] = w(arange(0,M+1),arange(0,C+1),hnvec)
    # Queue-length
    QN[arange[1,M+1],closedClasses] = multiply(matlib.repmat(XN(closedClasses),M,1),CN(arange(0,M+1),closedClasses))
    # Compute performance indexes at Nc for open classes
    for r in openClasses:
        # Throughput
        XN[r] = lambda_[r]
        for i in arange(0,M+1):
            # Queue-length
            QN[i,r] = 0
            for n in arange(0,sum(Nc)+1):
                QN[i,r] = QN[i,r] + lambda_[r] * D[i,r] * (n + 1) * EC(i,n + 1) * Pc(i,1 + n,hnvec)
            # Response time
            CN[i,r] = QN[i,r] / lambda_[r]
            # Utilization - the formula from Bruell-Balbo-Ashfari does not
            # match simulation, this appears to be simly lambda_r*D_{ir}
            UN[i,r] = 0
            for n in arange(0,sum(Nc)+1):
                UN[i,r] = UN[i,r] + lambda_[r] * Eprime(i,1 + n + 1) / E(i,1 + n + 1) * Pc(i,1 + n,hnvec)

    return XN,QN,UN,CN,lGN