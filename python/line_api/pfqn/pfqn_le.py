from numpy import *
from ..util import *
from ..util.allbut import allbut


def pfqn_le(L = None,N = None,Z = None):
    # [GN,LGN]=PFQN_LE(L,N,Z)

    # PFQN_LE Asymptotic solution of closed product-form queueing networks by
    # logistic expansion

    # [Gn,lGn]=pfqn_le(L,N,Z)
    # Input:
    # L : MxR demand matrix. L[i,r] is the demand of class-r at queue i
    # N : 1xR population vector. N[r] is the number of jobs in class r
    # Z : 1xR think time vector. Z[r] is the total think time of class r

    # Output:
    # Gn : estimated normalizing constat
    # lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
    #      will be correctly estimated.

    # Reference:
    # G. Casale. Accelerating performance inference over closed systems by
    # asymptotic methods. ACM SIGMETRICS 2017.
    # Available at: http://dl.acm.org/citation.cfm?id=3084445

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    if len(L)==0 or len(N)==0 or sum(N) == 0 or sum(L) < 0.0001:
        lGn = - sum(factln(N)) + sum(multiply(N,log(sum(Z, 0))))
        Gn = exp(lGn)
    else:
        if len(varargin) < 3:
            umax = pfqn_le_fpi(L,N)
            A = pfqn_le_hessian(L,N,transpose(umax))
            S = 0
            for r in arange(0,R+1):
                S = S + N[r] * log(transpose(umax) * L[:,r])
                lGn = multinomialln(array([N,M - 1])) + factln(M - 1) + (M - 1) * log(sqrt(2 * pi)) - log(sqrt(det(A))) + sum(log(umax)) + S
                Gn = exp(lGn)
            else:
                umax,vmax = pfqn_le_fpiZ(L,N,Z)
                A = pfqn_le_hessianZ(L,N,Z,transpose(umax),vmax)
                S = 0
                for r in arange(0,R+1):
                    S = S + N[r] * log(Z[r] + vmax * transpose(umax) * L[:,r])
                    lGn = - sum(factln(N)) - vmax + M * log(vmax) + M * log(sqrt(2 * pi)) - log(sqrt(det(A))) + sum(log(umax)) + S
                    Gn = exp(lGn)

            return Gn,lGn




def pfqn_le_fpi(L = None,N = None):
    # [U,D]=PFQN_LE_FPI(L,N)

    # find location of mode of gaussian
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    u = ones((M,1)) / M
    u_1 = Inf * u
    d = []
    while linalg.norm(u - u_1,1) > 1e-10:

        u_1 = u
        for i in arange(0,M+1):
            u[i] = 1 / (sum(N) + M)
            for r in arange(0,R+1):
                u[i] = u[i] + N[r] / (sum(N) + M) * L[i,r] * u_1[i] / (transpose(u_1) * L[:,r])
            d[end() + 1,:] = transpose(abs(u - u_1))


        return u,d


    def pfqn_le_fpiZ(L = None,N = None,Z = None):
        # [U,V,D]=PFQN_LE_FPIZ(L,N,Z)

        # find location of mode of gaussian
            if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
        eta = sum(N) + M
        u = ones((M,1)) / M
        v = eta + 1
        u_1 = Inf * u
        v_1 = Inf * v

        d = []
        while norm(u - u_1,1) > 1e-10:

            u_1 = u
            v_1 = v
            for i in arange(0,M+1):
                u[i] = 1 / eta
                for r in arange(0,R+1):
                    u[i] = u[i] + (N[r] / eta) * (Z[r] + v * L[i,r]) * u_1[i] / (Z[r] + v * transpose(u_1) * L[:,r])
                for r in arange(0,R+1):
                    xi[r] = N[r] / (Z[r] + v * transpose(u_1) * L[:,r])
                    v = eta + 1
                    for r in arange(0,R+1):
                        v = v - xi[r] * Z[r]
                    d[end() + 1,:] = transpose(abs(u - u_1)) + abs(v - v_1)


                return u,v,d


            def pfqn_le_hessian(L = None,N = None,u0 = None):
                # HU=PFQN_LE_HESSIAN(L,N,U0)

                # find hessian of gaussian
                    if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
                Ntot = sum(N)
                hu = zeros((M - 1,M - 1))
                for i in arange(0,(M - 1)+1):
                    for j in arange(0,(M - 1)+1):
                        if i != j:
                            hu[i,j] = - (Ntot + M) * u0[i] * u0[j]
                            for r in arange(0,R+1):
                                hu[i,j] = hu[i,j] + N[r] * L[i,r] * L[j,r] * (u0[i] * u0[j]) / (u0 * L[:,r]) ** 2
                        else:
                            hu[i,j] = (Ntot + M) * u0[i] * sum(allbut(u0,i))
                            for r in arange(0,R+1):
                                hu[i,j] = hu[i,j] - N[r] * L[i,r] * u0[i] * (allbut(u0,i) * L(allbut(arange(0,M+1),i),r)) / (u0 * L[:,r]) ** 2

                return hu


            def pfqn_le_hessianZ(L = None,N = None,Z = None,u = None,v = None):
                # A=PFQN_LE_HESSIANZ(L,N,Z,U,V)

                # find hessian of gaussian
                K,R = L.shape
                Ntot = sum(N)
                A = zeros((K,K))
                csi = zeros((1,R))
                for r in arange(0,R+1):
                    csi[r] = N[r] / (Z[r] + v * u * L[:,r])

                    Lhat = zeros([k,r])
                    for k in arange(0,K+1):
                        for r in arange(0,R+1):
                            Lhat[k,r] = Z[r] + v * L[k,r]

                    eta = Ntot + K
                    for i in arange(0,K+1):
                        for j in arange(0,K+1):
                            if i != j:
                                A[i,j] = - eta * u[i] * u[j]
                                for r in arange(0,R+1):
                                    A[i,j] = A[i,j] + csi[r] ** 2 * Lhat[i,r] * Lhat[j,r] * (u[i] * u[j]) / N[r]

                    for i in arange(0,K+1):
                        A[i,i] = - sum(allbut(A[i,:],i))

                        A = A[arange(0,(K - 1)+1),arange(0,(K - 1)+1)]
                        A[K,K] = 1
                        for r in arange(0,R+1):
                            A[K,K] = A(K,K) - (csi[r] ** 2 / N[r]) * Z[r] * u * L[:,r]

                            A[K,K] = v * A(K,K)
                            for i in arange(0,(K - 1)+1):
                                A[i,K] = 0
                                for r in arange(0,R+1):
                                    A[i,K] = A(i,K) + v * u[i] * ((csi[r] ** 2 / N[r]) * Lhat[i,r] * (u * L[:,r]) - csi[r] * L[i,r])
                                    A[K,i] = A(i,K)

                                return A


