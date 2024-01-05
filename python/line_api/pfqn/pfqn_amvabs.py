from numpy import *
import numpy.matlib as matlib

def pfqn_amvabs(L = None,N = None,Z = None,tol = None,maxiter = None,QN = None,weight = None):
    # [XN,QN,UN]=PFQN_AMVABS(L,N,Z,TOL,MAXITER,QN,WEIGHT)

    if Z is None:
        Z = 0 * N

    if tol is None:
        tol = 1e-06

    if maxiter is None:
        maxiter = 1000

    if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else:
        M = 0
        R = N.shape[1]

    if Z is None:
        Z = zeros((R))

    if ndim(Z) == 2:
        Z = sum(Z, 0)

    CN = zeros((M,R))
    if QN is None:
        QN = matlib.repmat(N,M,1) / M
    else:
        QN = QN + eps

    if weight is None:
        weight = ones((M,R))

    XN = zeros((1,R))
    UN = zeros((M,R))
    relprio = ones((M,R))
    for it in arange(0,maxiter+1):
        QN_1 = QN
        for i in arange(0,M):
            for r in arange(0,R):
                relprio[i,r] = (QN[i,r] * weight[i,r])
        for i in arange(0,M):
            for r in arange(0,R):
                CN[i,r] = L[i,r]
                for s in arange(0,R):
                    if s != r:
                        CN[i,r] = CN[i,r] + L[i,r] * QN[i,s] * relprio[i,s] / relprio[i,r]
                    else:
                        CN[i,r] = CN[i,r] + L[i,r] * QN[i,r] * (N[r] - 1) / N[r] * relprio[i,s] / relprio[i,r]
                XN[0,r] = N[r] / (Z[r] + sum(CN[:,r]))
        for i in arange(0,M):
            for r in arange(0,R):
                QN[i,r] = XN[0,r] * CN[i,r]
        for i in arange(0,M):
            for r in arange(0,R):
                UN[i,r] = XN[0,r] * L[i,r]
        if amax(abs(1 - QN / QN_1)) < tol:
            break

    return XN,QN,UN