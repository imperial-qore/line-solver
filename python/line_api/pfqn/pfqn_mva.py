from numpy import *


def pfqn_mva(L = None,N = None,Z = None,mi = None):
    # [XN,QN,UN,CN,LGN] = PFQN_MVA(L,N,Z,MI)

    if isoctave:
        #warning off;
        pass

    XN = []
    QN = []
    UN = []
    CN = []
    lGN = 0
    InfServ = 1
    if len(varargin) == 2:
        InfServ = 0

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]

    N = transpose(N)
    if len(varargin) < 4:
        mi = ones((1,M))

    if len(varargin) < 3 or len(Z)==0:
        Z = zeros((1,R))

    if (not any(N) ):
        #line_warning(mfilename,'closed populations are empty');
        return XN,QN,UN,CN,lGN

    NR = len(N)
    if (R != NR):
        line_error(mfilename,'demand matrix and population vector have different number of classes')

    XN = zeros((1,R))
    QN = zeros((M,R))
    CN = zeros((M,R))
    if InfServ == 1:
        Z = transpose(Z)
    else:
        Z = zeros((1,R))

    prods = zeros((1,R - 1))

    for w in arange(0,R - 1+1):
        prods[1,w] = prod(ones((1,R - (w + 1) + 1)) + N(1,arange(w + 1,R+1)))

    firstnonempty = R
    while (N(firstnonempty) == 0):

        firstnonempty = firstnonempty - 1


    totpop = prod(N + 1)
    ctr = totpop
    Q = zeros((totpop,M))
    currentpop = 2
    n = zeros((1,R))
    n[1,firstnonempty] = 1
    while ctr:

        s = 1
        while s <= R:

            pos_n_1s = 0
            if N[s] > 0:
                n[s] = N[s] - 1
                pos_n_1s = N[r]
                w = 1
                while w <= R - 1:

                    pos_n_1s = pos_n_1s + n(w) * prods(w)
                    w = w + 1

                n[s] = N[s] + 1
            CNtot = 0
            i = 1
            while i <= M:

                Lis = L(i,s)
                CN[i,s] = Lis * (mi[i] + Q(1 + pos_n_1s,i))
                CNtot = CNtot + CN(i,s)
                i = i + 1

            XN[s] = N[s] / (Z[s] + CNtot)
            i = 1
            while i <= M:

                QN[i,s] = XN[s] * CN(i,s)
                Q[currentpop,i] = Q(currentpop,i) + QN(i,s)
                i = i + 1

            s = s + 1

        s = R
        while s > 0 and (n(1,s) == N[s]) or s > firstnonempty:

            s = s - 1

        # now compute the normalizing constant
        last_nnz = find(n > 0,1,'last')
        if sum(n(arange(0,last_nnz - 1+1))) == sum(N(arange(0,last_nnz - 1+1))) and sum(n(arange((last_nnz + 1),R+1))) == 0:
            logX = log(XN(last_nnz))
            if not len(logX)==0 :
                lGN = lGN - logX
        if s == 0:
            break
        n[s] = N[s] + 1
        s = s + 1
        while s <= R:

            n[s] = 0
            s = s + 1

        ctr = ctr - 1
        currentpop = currentpop + 1


    for m in arange(0,M+1):
        for r in arange(0,R+1):
            UN[m,r] = XN[r] * L(m,r)

    return XN,QN,UN,CN,lGN