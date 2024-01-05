from numpy import *

def pfqn_fnc(alpha = None,c = None):
    # generate rates for functional server f(n)=n+c
    M = alpha.shape[0]
    if c is None:
        c = zeros((1,M))

    mu = pfqn_fnc(alpha,c)
    if not all(isfinite(mu)) :
        c = - 0.5 * ones((1,M))
        mu = pfqn_fnc(alpha,c)
    dt = 0
    it = 0
    while not all(isfinite(mu)) :
        it = it + 1
        dt = dt + 0.05
        c = - 1 / 2 + dt
        mu = pfqn_fnc(alpha,c)
        if c >= 2:
            break

    return mu,c

    N = len(alpha[1,:])
    mu = zeros((M,N))
    for i in arange(0,M+1):
        mu[i,1] = alpha(i,1) / (1 + c[i])
        for n in arange(2,N+1):
            rho = 0
            for k in arange(0,(n - 1)+1):
                rho = rho + (prod(alpha[i,arange((n - k + 1),n+1])) - prod(alpha[i,arange((n - k),(n - 1)+1)])) / prod(mu[i,arange(0,k+1)])
            mu[i,n] = (prod(alpha[i,arange(0,n+1)]) / prod(mu[i,arange(0,(n - 1)+1)]))
            mu[i,n] = mu[i,n] / (1 - rho)

    mu[isnan(mu)] = Inf
    mu[abs(mu) > 1000000000000000.0] = Inf
    for i in arange(0,M+1):
        if any(isinf(mu[i,:])):
            s = amin(nonzero(isinf(mu[i,:])))
            mu[i,arange(s, mu.shape[1])] = Inf

#mu(mu==0) = Inf;
return mu,c