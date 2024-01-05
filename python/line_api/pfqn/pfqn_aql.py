from numpy import *
from ..util import *

def pfqn_aql(L = None,N = None,Z = None,TOL = None,QN0 = None):

    MAXITER = 1000.0

    K,C = L.shape
    if Z is None:
        Z = zeros((1,C))

    if TOL is None:
        TOL = 1e-07

    Q = []
    for k in arange(0,C):
        Qt = ones(K)
        Q.append(Qt)

    R = []
    for k in arange(0,C):
        Rt = ones((K,C))
        R.append(Rt)

    X = []
    for k in arange(0,C):
        Xt = ones(C)
        X.append(Xt)

    gamma = zeros((K,C))

    if QN0 is None:
        for t in arange(0,C):
            n = oner(N,t)
            for k in arange(0,K):
                Q[t][k] = sum(n) / K
    else:
        for t in arange(0,C):
            n = oner(N,t)
            for k in arange(0,K):
                Q[t][k] = QN0[k]

    it = 0
    while 1:
        Q_olditer = Q
        it = it + 1
        for t in arange(0,C):
            n = oner(N,t)
            for k in arange(0,K):
                for s in arange(0,C):
                    R[t][k,s] = L[k,s] * (1 + (sum(n) - 1) * (Q[t][k] / sum(n) - gamma[k,s]))
            for s in arange(0,C):
                X[t][s] = n[s] / (Z[s] + sum(R[t][:,s]))
            for k in arange(0,K):
                Q[t][k] = transpose(X[t]) @ transpose(R[t][k,:])
        for k in arange(0,K):
            for s in arange(0,C):
                gamma[k,s] = (Q[0][k] / sum(N)) - (Q[s][k] / (sum(N) - 1))
        if amax(abs((Q_olditer[0] - Q[0]) / Q[0])) < TOL or it == MAXITER:
            numIters = it
            break

    XN = X[0]
    RN = R[0]
    UN = zeros((K,C))
    QN = zeros((K,C))
    AN = zeros((K,C))
    for k in arange(0,K):
        for s in arange(0,C):
            UN[k,s] = XN[s] * L[k,s]
            QN[k,s] = UN[k,s] * (1 + Q[s][k])
            AN[k,s] = Q[s][k]

    return XN,QN,UN,RN,AN,numIters
