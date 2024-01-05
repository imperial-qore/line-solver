from numpy import *
from ..util import *
def pfqn_procomom2(L = None,N = None,Z = None,mu = None,m = None):
    # Marginal state probabilities for the queue in a model consisting of a
    # queueing station and a delay station only.

    if len(varargin) < 4 or len(mu)==0:
        mu = ones((m,sum(N) + 1))
    else:
        mu = array([1,transpose(mu)])

    if len(varargin) < 5:
        m = 1

    __,R = L.shape
    # compute solution for [1,0,0,...,0]
    p0 = zeros((sum(N) + 1,1))
    p0[end()] = 1
    # compute the rest
    tic
    for r in arange(0,R+1):
        # generate F2r matrix
        T[r] = sparse(1 + sum(N),1 + sum(N))
        for n in arange(sum(N),1+- 1,- 1):
            row = sum(N) - n + 1
            T[r][row,row] = Z[r]
            T[r][row,row + 1] = (n + m - 1) * L[r] / mu(1 + n)
        T[r][sum[N] + 1,sum[N] + 1] = Z[r]

    F = eye(sum(N) + 1)
    B = eye(sum(N) + 1)
    for r in arange(0,R+1):
        F = F * T[r] ** N[r] / factorial(N[r])
        B = B * T[r]

    pk = transpose((F * p0))
    G = sum(pk)
    if any(not isfinite(pk(1,1)) ) or not isfinite(G) :
        # todo
        lG = logsumexp(log(pk(1,:)))
        else:
        if not isfinite(G) :
            lG = logsumexp(log(pk(1,:)))
            else:
            lG = log(G)

    pk = pk / G
    pk = pk(arange(end(),1+- 1,- 1))
    ## test
    Q = 0
    # V=0;
    for n in arange(0,sum(N)+1):
        psingle[n + 1] = pk(n + 1)

    psingle = psingle / sum(psingle)
    for n in arange(0,sum(N)+1):
        Q = Q + n * psingle(n + 1)
        #     V=V + (n^2-n) * psingle(n+1);

    # psingle
    QN = double(Q)
    XNMVA,QNMVA,__,__,__,__,pik = pfqn_mvald(matlib.repmat(L,m,1),N,Z,matlib.repmat(mu(:,arange(2,-1)),m,1))
    # VNMVA=0;
    # for n=1:sum(N)
    #     VNMVA=VNMVA+(n^2-n)*pik(1,:);
    # end
    # QNMVA = sum(QNMVA,2)
    # QNTOTMVA=sum(QNMVA)
    # [V,VNMVA]
    return pk,lG,G,T,F,B