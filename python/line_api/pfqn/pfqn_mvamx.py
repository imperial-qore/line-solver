from numpy import *

def pfqn_mvamx(lambda_ = None,D = None,N = None,Z = None,mi = None):
    # [XN,QN,UN,CN,LGN] = PFQN_MVAMX(LAMBDA,D,N,Z, MI)

    if any(N(lambda_ > 0) > logical_and(0,isfinite(N(lambda_ > 0)))):
        line_error(mfilename,'Arrival rate cannot be specified on closed classes.')

    M,R = D.shape
    if len(varargin) < 5:
        mi = ones((M,1))

    openClasses = find(isinf(N))
    closedClasses = setdiff1d(arange(0,len(N)+1),openClasses)
    XN = zeros((1,R))
    UN = zeros((M,R))
    CN = zeros((M,R))
    QN = zeros((M,R))
    for r in openClasses:
        for i in arange(0,M+1):
            UN[i,r] = lambda_[r] * D[i,r]
        XN[r] = lambda_[r]

    UNt = sum(UN, 1)
    if len(Z)==0:
        Z = zeros((1,R))

    Dc = D(:,closedClasses) / (1 - matlib.repmat(UNt,1,len((closedClasses))))
    XNc,QNc,__,CNc,lGN = pfqn_mva(Dc,N(closedClasses),Z(closedClasses),mi)
    XN[closedClasses] = XNc
    QN[:,closedClasses] = QNc
    CN[:,closedClasses] = CNc
    for i in arange(0,M+1):
        for r in closedClasses:
            UN[i,r] = XN[r] * D[i,r]

    for i in arange(0,M+1):
        for r in openClasses:
            if len(QNc)==0:
                CN[i,r] = D[i,r] / (1 - UNt[i])
            else:
                CN[i,r] = D[i,r] * (1 + sum(QNc[i,:])) / (1 - UNt[i])
            QN[i,r] = CN[i,r] * XN[r]

    return XN,QN,UN,CN,lGN