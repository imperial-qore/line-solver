from numpy import *
from ..util import *

def pfqn_mvaldmx_ec(lambda_ = None,D = None,mu = None):
    # [EC,E,EPRIME,LO] = PFQN_MVALDMX_EC(LAMBDA,D,MU)
    # Compute the effective capacity terms in MVALDMX
    # Think times are not handled since this assumes limited load-dependence
    M,R = mu.shape
    #Nt = sum(N(isfinite(N)));
    Lo = zeros((M,1))
    for i in arange(0,M+1):
        Lo[i] = lambda_ * transpose(D[i,:])

        b = zeros((M,1))

        for i in arange(0,M+1):
            b[i] = find(mu[i,:] == mu(i,end()),1)

            Nt = mu.shape[1]

            mu[:,arange[end() + 1,end() + 1 + amax[b]+1]] = matlib.repmat(mu(:,end()),1,1 + amax(b))
            C = 1.0 / mu
            EC = zeros((M,Nt))
            #Ever = zeros(M,1+Nt);
            E = zeros((M,1 + Nt))
            Eprime = zeros((M,1 + Nt))
            for i in arange(0,M+1):
                E1 = zeros((1 + Nt,1 + Nt))
                E2 = zeros((1 + Nt,1 + Nt))
                E3 = zeros((1 + Nt,1 + Nt))
                F2 = zeros((1 + Nt,1 + b[i] - 2))
                F3 = zeros((1 + Nt,1 + b[i] - 2))
                E2prime = zeros((1 + Nt,1 + Nt))
                F2prime = zeros((1 + Nt,1 + b[i] - 2))
                for n in arange(0,Nt+1):
                    if n >= b[i]:
                        E[i,1 + n] = 1 / (1 - Lo[i] * C(i,b[i])) ** (n + 1)
                        #            Ever(i,1+n) = 1 / (1-Lo[i]*C(i,b[i]))^(n+1);
                        Eprime[i,1 + n] = C(i,b[i]) * E(i,1 + n)
                    else:
                        ## compute E1
                        if n == 0:
                            E1[1 + n] = 1 / (1 - Lo[i] * C(i,b[i]))
                            for j in arange(0,(b[i] - 1)+1):
                                E1[1 + n] = E1(1 + n) * C[i,j] / C(i,b[i])
                        else:
                            E1[1 + n] = 1 / (1 - Lo[i] * C(i,b[i])) * C(i,b[i]) / C(i,n) * E1(1 + (n - 1))
                        ## compute F2
                        for n0 in arange(0,(b[i] - 2)+1):
                            if n0 == 0:
                                F2[1 + n,1 + n0] = 1
                            else:
                                F2[1 + n,1 + n0] = (n + n0) / n0 * Lo[i] * C(i,n + n0) * F2(1 + n,1 + (n0 - 1))
                        ## compute E2
                        E2[1 + n] = sum(F2(1 + n,1 + (arange(0,b[i] - 2+1))))
                        ## compute F3
                        for n0 in arange(0,(b[i] - 2)+1):
                            if n == 0 and n0 == 0:
                                F3[1 + n,1 + n0] = 1
                                for j in arange(0,(b[i] - 1)+1):
                                    F3[1 + n,1 + n0] = F3(1 + n,1 + n0) * C[i,j] / C(i,b[i])
                            else:
                                if n > 0 and n0 == 0:
                                    F3[1 + n,1 + n0] = C(i,b[i]) / C(i,n) * F3(1 + (n - 1),1 + 0)
                                else:
                                    F3[1 + n,1 + n0] = (n + n0) / n0 * Lo[i] * C(i,b[i]) * F3(1 + n,1 + (n0 - 1))
                        ## compute E3
                        E3[1 + n] = sum(F3(1 + n,1 + (arange(0,b[i] - 2+1))))
                        ## compute F2prime
                        for n0 in arange(0,(b[i] - 2)+1):
                            if n0 == 0:
                                F2prime[1 + n,1 + n0] = C(i,n + 1)
                            else:
                                F2prime[1 + n,1 + n0] = (n + n0) / n0 * Lo[i] * C(i,n + n0 + 1) * F2prime(1 + n,1 + (n0 - 1))
                        ## compute E2prime
                        E2prime[1 + n] = sum(F2prime(1 + n,1 + (arange(0,(b[i] - 2)+1))))
                        # finally, compute E, Eprime, and EC
                        E[i,1 + n] = E1(1 + n) + E2(1 + n) - E3(1 + n)
                        if n < b[i] - 1:
                            Eprime[i,1 + n] = C(i,b[i]) * E1(1 + n) + E2prime(1 + n) - C(i,b[i]) * E3(1 + n)
                        else:
                            Eprime[i,1 + n] = C(i,b[i]) * E(i,1 + n)
                        ## verification of E
                        #            Ever(i,1+n) = C(i,b[i])^(n+1-b[i]) / (1-Lo[i]*C(i,b[i]))^(n+1) * prod(C(i,(n+1):(b[i]-1)));
                        #            for n0 = 0:(b[i]-2)
                        #                Ever(i,1+n) = Ever(i,1+n) + nchoosek(n+n0,n0) * Lo[i]^n0 * (prod(C(i,(n+1):(n+n0)))-C(i,b[i])^(n0+n+1-b[i])*prod(C(i,(n+1):(b[i]-1))));
                        #            end
                        #         ## verification2 of E
                        #         Ever(i,1+n) = 0;
                        #         for n0=0:1000
                        #             if n+n0+1>b[i]
                        #                 C(i,n+n0+1) = 1/mu(i,b[i]);
                        #             end
                        #             Ever(i,1+n) = Ever(i,1+n)  +nchoosek(n+n0,n0) * Lo[i]^n0 * prod(C(i,(n+1):(n+n0)));
                        #         end
                        #         ## verification of Eprime
                        #         Eprimever(i,1+n) = 0;
                        #         for n0=0:1000
                        #             if n+n0+1>b[i]
                        #                 C(i,n+n0+1) = 1/mu(i,b[i]);
                        #             end
                        #             Eprimever(i,1+n) = Eprimever(i,1+n)  +nchoosek(n+n0,n0) * Lo[i]^n0 * prod(C(i,(n+1):(n+n0+1)));
                        #         end
                        #    Eprime = Eprimever;
                        # EC not defined for n=0
                for n in arange(0,Nt+1):
                    #if n>=b[i]
                    #    EC(i,n) = C(i,b[i]) / (1-Lo[i]*C(i,b[i]));
                    #elseif n>0
                    EC[i,n] = C(i,n) * E(i,1 + n) / E(i,1 + (n - 1))
                    #end
                    #EC
                    #E
                    #Ever
                    #Eprime
                    #Eprimever
                    #Eprime = Eprimever;
            return EC,E,Eprime,Lo