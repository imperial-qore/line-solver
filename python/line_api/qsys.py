import numpy as np
    
def qsys_gig1_approx_allencunneen(lambda_ = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_ALLENCUNNEEN(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    W = (rho / (1 - rho)) / mu * ((cs ** 2 + ca ** 2) / 2) + 1 / mu
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gig1_approx_gelenbe(lambda_ = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_GELENBE(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    W = (rho * ca ** 2 + cs ** 2) / 2 / (1 - rho) / lambda_
    return W
    
def qsys_gig1_approx_heyman(lambda_ = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_HEYMAN(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    W = rho / (1 - rho) / mu * (ca ** 2 + cs ** 2) / 2 + 1 / mu
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gig1_approx_kimura(sigma = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_KIMURA(SIGMA,MU,CA,CS)
    
    W = sigma * (ca ** 2 + cs ** 2) / mu / (1 - sigma) / (1 + ca ** 2)
    return W
    
def qsys_gig1_approx_klb(lambda_ = None,mu = None,ca = None,cs = None): 
    # [W,rhohat]=QSYS_GIG1_APPROX_KLB(LAMBDA,MU,CA,CS)
    
    # kramer-langenbach-belz formula
    rho = lambda_ / mu
    if ca <= 1:
        g = np.exp(- 2 * (1 - rho) * (1 - ca ** 2) ** 2 / (3 * rho * (ca ** 2 + cs ** 2)))
    else:
        g = np.exp(- (1 - rho) * (ca ** 2 - 1) / (ca ** 2 + 4 * cs ** 2))
    
    W = 1 / mu * ((rho / (1 - rho)) * ((cs ** 2 + ca ** 2) / 2) * g + 1)
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gig1_approx_kobayashi(lambda_ = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_KOBAYASHI(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    rhohat = np.exp(- 2 * (1 - rho) / (rho * (ca ** 2 + cs ** 2 / rho)))
    W = rhohat / (1 - rhohat) / lambda_
    return W,rhohat
    
def qsys_gig1_approx_marchal(lambda_ = None,mu = None,ca = None,cs = None): 
    # W=QSYS_GIG1_APPROX_MARCHAL(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    Wmm1 = rho / (1 - rho)
    W = Wmm1 * (1 + cs ** 2) / 2 / mu * (ca + rho ** 2 * cs ** 2) / (1 + rho ** 2 * cs ** 2) + 1 / mu
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gig1_approx_myskja(lambda_ = None,mu = None,ca = None,cs = None,q0 = None,qa = None): 
    # W=QSYS_GIG1_APPROX_MYSKJA(LAMBDA,MU,CA,CS,Q0,QA)
    
    # qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
    # q0 = lowest value of the relative third moment for a given mean and SCV
    rho = lambda_ / mu
    W = rho / (2 * mu * (1 - rho)) * ((1 + cs ** 2) + (q0 / q) ** (1 / rho - rho) * (1 / rho) * (ca ** 2 - 1))
    return W
    
def qsys_gig1_approx_myskja2(lambda_ = None,mu = None,ca = None,cs = None,q0 = None,qa = None): 
    # W=QSYS_GIG1_APPROX_MYSKJA2(LAMBDA,MU,CA,CS,Q0,QA)
    
    # qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
    # q0 = lowest value of the relative third moment for a given mean and SCV
    ra = (1 + ca ** 2) / 2
    rs = (1 + cs ** 2) / 2
    rho = lambda_ / mu
    theta = (rho * (qa - ra) - (qa - ra ** 2)) / (2 * rho * (ra - 1))
    d = (1 + 1 / ra) * (1 - rs) * (1 - (q0 / qa) ** 3) * (1 - rho ** 3)
    D = (rs - theta) ** 2 + (2 * rs - 1 + d) * (ra - 1)
    W = (rho / (1 - rho)) / lambda_ * (rs + (1 / rho) * (np.sqrt(D) - (rs - theta)))
    
    return W
    
def qsys_gig1_ubnd_kingman(lambda_ = None,mu = None,ca = None,cs = None): 
    # WUB=QSYS_GIG1_UBND_KINGMAN(LAMBDA,MU,CA,CS)
    
    rho = lambda_ / mu
    W = rho / (1 - rho) * (ca ** 2 + cs ** 2) / 2 * (1 / mu) + (1 / mu)
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gigk_approx(lambda_ = None,mu = None,ca = None,cs = None,k = None): 
    # W=QSYS_GIGK_APPROX(LAMBDA,MU,CA,CS,K)
    rho = lambda_ / (mu * k)
    if rho > 0.7:
        alpha = (rho ** k + rho) / 2
    else:
        alpha = rho ** ((k + 1) / 2)
    
    W = (alpha / mu) * (1 / (1 - rho)) * (ca ** 2 + cs ** 2) / (2 * k) + 1 / mu
    
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gigk_approx_kingman(lambda_ = None,mu = None,ca = None,cs = None,k = None): 
    # W=QSYS_GIG1_UBND_KINGMAN(LAMBDA,MU,CA,CS,K)
    W = (ca ** 2 + cs ** 2) / 2 * (qsys_mmk(lambda_,mu,k) - 1 / mu) + 1 / mu
    rhohat = W * lambda_ / (1 + W * lambda_)
    
    return W,rhohat
    
def qsys_gm1(sigma = None,mu = None): 
    # W=QSYS_GM1(SIGMA,MU)
    
    # sigma = Load at arrival instants (Laplace transform of the inter-arrival times)
    W = 1 / (1 - sigma) / mu
    
    return W
    
def qsys_mg1(lambda_ = None,mu = None,cs = None): 
    # W=QSYS_MG1(LAMBDA,MU,CS)
    
    rho = lambda_ / mu
    Q = rho + rho ** 2 / (2 * (1 - rho)) + lambda_ ** 2 * cs ** 2 / mu ** 2 / (2 * (1 - rho))
    W = Q / lambda_
    rhohat = Q / (1 + Q)
    
    return W,rhohat
    
def qsys_mm1(lambda_ = None,mu = None): 
    # W=QSYS_MM1(LAMBDA,MU)
    rho = lambda_ / mu
    W = rho / (1 - rho) / lambda_
    
    return W,rho
    
def qsys_mmk(lambda_ = None,mu = None,k = None): 
    # W=QSYS_MMK(LAMBDA,MU,K)
    rho = lambda_ / mu / k
    Q = rho / (1 - rho) * qsys_erlangc(k,rho) + k * rho
    W = Q / lambda_
    return W,rho
    
def qsys_erlangc(k = None,rho = None): 
    # Erlang-C formula
    # The probability that an arriving customer is forced to join the queue
    # (all servers are occupied)
    S = 0
    for j in np.arange(0,(k - 1)+1).reshape(-1):
        S = S + (k * rho) ** j / factorial[j]
    
    C = 1 / (1 + (1 - rho) * factorial[k] / (k * rho) ** k * S)
    return C
