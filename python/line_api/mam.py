import numpy as np
    
def qbd_R(B = None,L = None,F = None,iter_max = None): 
    # Successive substitutions method
    if len(varargin) < 4:
        iter_max = 100000
    
    Fil = F * inv(L)
    BiL = B * inv(L)
    R = - Fil
    Rprime = - Fil - R ** 2 * BiL
    for iter in np.arange(1,iter_max+1).reshape(-1):
        R = Rprime
        Rprime = - Fil - R ** 2 * BiL
        if norm(R - Rprime,1) <= 1e-12:
            break
    
    R = Rprime
    return R
    
def qbd_R_logred(B = None,L = None,F = None,iter_max = None): 
    # Logarithmic reduction method
    if len(varargin) < 4:
        iter_max = 100000
    
    iLF = - inv(L) * F
    iLB = - inv(L) * B
    T = iLF
    S = iLB
    for iter in np.arange(1,iter_max+1).reshape(-1):
        D = iLF * iLB + iLB * iLF
        iLF = inv(np.eye[r] - D) * iLF * iLF
        iLB = inv(np.eye[r] - D) * iLB * iLB
        S = S + T * iLB
        T = T * iLF
        if norm(np.ones((r,1)) - S * np.ones((r,1)),1) <= 1e-12:
            break
    
    U = L + F * S
    R = - F * inv(U)
    return R
