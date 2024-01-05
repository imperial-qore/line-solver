from numpy import *
from ..util import *

def pfqn_rd(L = None,N = None,Z = None,mu = None,options = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    lambda_ = zeros((1,R))
    if sum(N) < 0:
        lGN = - Inf
        return lGN,Cgamma

    if len(varargin) < 5:
        options = SolverNC.defaultOptions

    # L
    # N
    # Z
    # mu
    for i in arange(0,M+1):
        if all(mu[i,:] == mu(i,1)):
        L[i,:] = L[i,:] / mu(i,1)
        mu[i,:] = 1
        #isLI[i] = true;

if sum(N) == 0:
    lGN = 0
    return lGN,Cgamma

gamma = ones((M,sum(N)))
mu = mu(:,arange(0,sum(N)+1))
#mu(mu==0)=Inf;
mu[isnan[mu]] = Inf
s = zeros((M,1))
for i in arange(0,M+1):
    if isfinite(mu(i,end())):
        s[i] = amin(find(abs(mu[i,:] - mu(i,end())) < options.tol))
        if s[i] == 0:
            s[i] = sum(N)
    else:
        s[i] = sum(N)

isDelay = False(M,1)
isLI = False(M,1)
y = L
for i in arange(0,M+1):
    if isinf(mu(i,s[i])):
        lastfinite = amax(find(isfinite(mu[i,:])))
        s[i] = lastfinite
    y[i,:] = y[i,:] / mu(i,s[i])

for i in arange(0,M+1):
    gamma[i,:] = mu[i,:] / mu(i,s[i])
    if amax(abs(mu[i,:] - (arange(0,sum(N)+1)))) < options.tol:
    #isDelay[i] = true;
    pass

# eliminating the delays seems to produce problems
# Z = sum([Z; L(isDelay,:)],1);
# L(isDelay,:)=[];
# mu(isDelay,:)=[];
# gamma(isDelay,:)=[];
# y(isDelay,:)=[];
# isLI(isDelay) = [];
# M = M - sum(isDelay);

beta = ones((M,sum(N)))
for i in arange(0,M+1):
    beta[i,1] = scipy.special.gamma(i,1) / (1 - scipy.special.gamma(i,1))
    for j in arange(2,sum(N)+1):
        beta[i,j] = (1 - scipy.special.gamma(i,j - 1)) * (scipy.special.gamma[i,j] / (1 - scipy.special.gamma[i,j]))

beta[isnan[beta]] = Inf
if (all(beta == Inf)):
    options.method = 'default'
    lGN = pfqn_nc(lambda_,L,N,Z,options)
    return lGN,Cgamma
else:
    Cgamma = 0
    sld = s(s > 1)
    vmax = amin(sum(sld - 1),sum(N))
    Y = pfqn_mva(y,N,0 * N)
    rhoN = y * transpose(Y)
    for vtot in arange(0,vmax+1):
        lEN[vtot + 1] = pfqn_gldsingle(rhoN,vtot,beta)
    lEN = real(lEN)
    for vtot in arange(0,vmax+1):
        EN = exp(lEN(vtot + 1))
        Cgamma = Cgamma + ((sum(N) - amax(0,amax(vtot - 1))) / sum(N)) * EN
    options.method = 'default'
    lGN = pfqn_nc(lambda_,y,N,Z,options)
    lGN = lGN + log(Cgamma)

return lGN,Cgamma