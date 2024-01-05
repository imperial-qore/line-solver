from numpy import *
from ..util import *

def pfqn_schmidt(D = None,N = None,S = None,sched = None):
    # [XN,QN,UN,CN] = PFQN_SCHMIDT(D,N,S,SCHED)

    # utilization in general ld case does not work
    M,R = D.shape
    closedClasses = arange(0,R+1)
    XN = zeros((1,R))
    UN = zeros((M,R))
    CN = zeros((M,R))
    QN = zeros((M,R))
    C = len(closedClasses)

    Dc = D[:,closedClasses]
    Nc = N[closedClasses]
    prods = zeros((1,C))

    for r in arange(0,C+1):
        prods[r] = prod(Nc(arange(0,r - 1+1)) + 1)

    # Start at nc=(0,...,0)
    kvec = pprod(Nc)
    # Initialize L and Pc
    L = array([])

    Pc = array([])

    for i in arange(0,M+1):
        if SchedStrategy.ID_INF == sched[i]:
            L[i] = zeros((R,prod(1 + Nc)))
        else:
            if SchedStrategy.ID_PS == sched[i]:
                if S[i] == 1:
                    L[i] = zeros((R,prod(1 + Nc)))
                else:
                    Pc[i] = zeros((1 + sum(Nc),prod(1 + Nc)))
            else:
                if SchedStrategy.ID_FCFS == sched[i]:
                    if all(D[i,:] == D(i,1)):
                    if S[i] == 1:
                        L[i] = zeros((R,prod(1 + Nc)))
                    else:
                        Pc[i] = zeros((1 + sum(Nc),prod(1 + Nc)))
                else:
                    Pc[i] = zeros((prod(1 + Nc),prod(1 + Nc)))

x = zeros((C,prod(1 + Nc)))
w = zeros((M,C,prod(1 + Nc)))
for i in arange(0,M+1):
    Pc[i][1 + 0,hashpop[kvec,Nc,C,prods]] = 1.0

u = zeros((M,C))
# Population recursion
while kvec >= 0:

    hkvec = hashpop(kvec,Nc,C,prods)
    nc = sum(kvec)
    kprods = zeros((1,C))
    for r in arange(0,C+1):
        kprods[r] = prod(kvec(arange(0,r - 1+1)) + 1)
    for i in arange(0,M+1):
        for c in arange(0,C+1):
            hkvec_c = hashpop(oner(kvec,c),Nc,C,prods)
            # Compute mean residence times
            for n in arange(0,nc+1):
                if SchedStrategy.ID_INF == sched[i]:
                    w[i,c,hkvec] = D(i,c)
                else:
                    if SchedStrategy.ID_PS == sched[i]:
                        if S[i] == 1:
                            w[i,c,hkvec] = Dc(i,c) * (1 + L[i](1 + n,hkvec_c))
                        else:
                            w[i,c,hkvec] = (Dc(i,c) / S[i]) * (1 + L[i](c,hkvec_c))
                            for i in arange(0,S[i] - 2+1):
                                w[i,c,hkvec] = w(i,c,hkvec) + (S[i] - 1 - i) * Pc[i](1 + i,hkvec_c)
                    else:
                        if SchedStrategy.ID_FCFS == sched[i]:
                            if all(D[i,:] == D(i,1)):
                            if S[i] == 1:
                                w[i,c,hkvec] = Dc(i,c) * (1 + L[i](1 + n,hkvec_c))
                            else:
                                w[i,c,hkvec] = (Dc(i,c) / S[i]) * (1 + L[i](c,hkvec_c))
                                for i in arange(0,S[i] - 2+1):
                                    w[i,c,hkvec] = w(i,c,hkvec) + (S[i] - 1 - i) * Pc[i](1 + i,hkvec_c)
                        else:
                            nvec = pprod(kvec)
                            while nvec >= 0:

                                if sum(nvec) > 0:
                                    hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods)
                                    Bcn = D(i,c) + amax(0,sum(nvec) - S[i]) / (S[i] * (sum(nvec) - 1)) * (nvec * transpose(D[i,:]) - D(i,c))
                                    w[i,c,hnvec] = w(i,c,hnvec) + Bcn * Pc[i](hnvec_c,hkvec_c)
                                nvec = pprod(nvec,kvec)

# Compute tput
for c in arange(0,C+1):
    x[c,hkvec] = kvec(c) / sum(w(arange(0,M+1),c,hkvec))
for i in arange(0,M+1):
    for c in arange(0,C+1):
        L[i][c] = x(c,hkvec) * w(i,c,hkvec)
    if SchedStrategy.ID_PS == sched[i]:
        if S[i] > 1:
            for n in arange(0,amin(S[i],sum(kvec))+1):
                for c in arange(0,C+1):
                    hkvec_c = hashpop(oner(kvec,c),Nc,C,prods)
                    Pc[i][1 + n,hkvec] = Pc[i](1 + n,hkvec) + Dc(i,c) * (1 / n) * x(c,hkvec) * Pc[i](1 + (n - 1),hkvec_c)
            Pc[i][1 + 0,hkvec] = amax(eps,1 - sum(Pc(i,1 + (arange(0,amin(S[i],sum(kvec))+1)),hkvec)))
    else:
        if SchedStrategy.ID_FCFS == sched[i]:
            if all(D[i,:] == D(i,1)):
            if S[i] > 1:
                for n in arange(0,amin(S[i],sum(kvec))+1):
                    for c in arange(0,C+1):
                        hkvec_c = hashpop(oner(kvec,c),Nc,C,prods)
                        Pc[i][1 + n,hkvec] = Pc[i](1 + n,hkvec) + Dc(i,c) * (1 / n) * x(c,hkvec) * Pc[i](1 + (n - 1),hkvec_c)
                if sum(kvec) > 0:
                    Pc[i][1 + 0,hkvec] = amax(eps,1 - sum(Pc(i,1 + (arange(0,amin(S[i],sum(kvec))+1)),hkvec)))
        else:
            nvec = pprod(kvec)
            while nvec >= 0:

                hnvec = hashpop(nvec,kvec,C,kprods)
                if sum(nvec) > 0:
                    for c in arange(0,C+1):
                        hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods)
                        hkvec_c = hashpop(oner(kvec,c),Nc,C,prods)
                        Bcn = D(i,c) + amax(0,sum(nvec) - S[i]) / (S[i] * (sum(nvec) - 1)) * (nvec * transpose(D[i,:]) - D(i,c))
                        Pc[i][hnvec,hkvec] = Pc[i](hnvec,hkvec) + (1 / nvec(c)) * x(c,hkvec) * Bcn * Pc[i](hnvec_c,hkvec_c)
                nvec = pprod(nvec,kvec)

kvec = pprod(kvec,Nc)


return XN,QN,UN,CN
