from numpy import *

def pfqn_propfair(L = None,N = None,Z = None):
    # [G,LOG] = PFQN_PROPFAIR(L,N,Z)

    # Proportionally fair allocation

    # Estimate the normalizing constant using a convex optimization program
    # that is asymptotically exact in models with single-server PS queues only.
    # The underlying optimization program is convex.
    # The script implements a heuristic to estimate the solution in the
    # presence of delay stations.

    # Schweitzer, P. J. (1979). Approximate analysis of multiclass closed networks of
    # queues. In Proceedings of the International Conference on Stochastic Control
    # and Optimization. Free Univ., Amsterdam.

    # Walton, Proportional fairness and its relationship with multi-class
    # queueing networks, 2009.

        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    optimopt = optimoptions(fmincon,'MaxFunctionEvaluations',1000000.0,'Display','none')
    obj = lambda X = None: - sum(multiply(N,log(X + 1e-06)))
    x0 = N / (multiply(N,sum(L, 0)) + Z)
    Xopt = fmincon(lambda x = None: obj(x),x0,L,ones((M,1)),[],[],zeros((1,R)),[],[],optimopt)
    Xasy = multiply(Xopt,N) / (N + multiply(Z,Xopt))

    lG = - sum(multiply(N,log(Xasy)))
    G = exp(lG)
    return G,lG,Xasy