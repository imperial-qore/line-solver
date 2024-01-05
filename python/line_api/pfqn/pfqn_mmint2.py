from numpy import *
from ..util import *

def pfqn_mmint2(L = None,N = None,Z = None,method = None):
    # [G,LOGG] = PFQN_PNC2(L,N,Z)

    if len(varargin) < 4:
        method = 'default'

    nnzClasses = find(N)
    # repairmen integration
    order = 12
    # below we use a variable substitution u->u^2 as it is numerically better
    if 'quadratic' == method:
        f = lambda u = None: transpose((multiply(multiply(2 * transpose(u),exp(- transpose((u ** 2)))),prod((Z(nnzClasses) + multiply(L(nnzClasses),matlib.repmat(u ** 2,1,len(nnzClasses)))) ** N(nnzClasses), 1))))
    else:
        if 'default' == method:
            f = lambda u = None: transpose((multiply(exp(- transpose(u)),prod((Z(nnzClasses) + multiply(L(nnzClasses),matlib.repmat(u,1,len(nnzClasses)))) ** N(nnzClasses), 1))))

    p = 1 - 10 ** - order
    exp1prctile = - log(1 - p) / 1

    w = warning
    warnings.warn('off')
    lG = log(integral(f,0,exp1prctile,'AbsTol',10 ** - order)) - sum(factln(N))
    G = exp(lG)
    warnings.warn(w)
    return G,lG