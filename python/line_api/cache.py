import numpy as np
import numpy.matlib
import scipy.special
    
def cache_erec(gamma = None,m = None): 
    if len(varargin) < 3:
        R = np.ones((len(m),len(gamma)))
    
    E = sub_cache_erec(gamma,m,len(gamma))
    return E
    
    
def sub_cache_erec(gamma = None,m = None,k = None): 
    h = len(m)
    if sum(m) == 0:
        E = 1
        return E
    
    if sum(m) > k or np.amin(m) < 0:
        E = 0
        return E
    
    if k == 1 and sum(m) == 1:
        j = find(m)
        E = scipy.special.gamma(1,j)
        return E
    
    E = sub_cache_erec(gamma,m,k - 1)
    for j in np.arange(1,h+1).reshape(-1):
        if m[j] > 0:
            E = E + scipy.special.gamma[k,j] * m[j] * sub_cache_erec(gamma,oner(m,j),k - 1)
    
    return E
    
    
def cache_gamma(lambda_ = None,R = None): 
    u = lambda_.shape[1-1]
    
    n = lambda_.shape[2-1]
    
    h = lambda_.shape[3-1] - 1
    
    gamma = np.zeros((n,h))
    for i in np.arange(1,n+1).reshape(-1):
        for j in np.arange(1,h+1).reshape(-1):
            ##
# compute gamma[i,j]
            G = digraph(R[1,i])
            Pj = G.shortestpath(1,j)
            gamma[i,j] = sum(lambda_(:,i,1 + 0))
            for li in np.arange(2,len(Pj)+1).reshape(-1):
                y = 0
                l_1 = Pj(li - 1)
                l = Pj(li)
                for v in np.arange(1,u+1).reshape(-1):
                    y = y + lambda_(v,i,1 + l_1) * R[v,i](l_1,l)
                gamma[i,j] = scipy.special.gamma[i,j] * y
    
    return gamma,u,n,h
    
def cache_gamma_lp(lambda_ = None,R = None): 
    u = lambda_.shape[1-1]
    
    n = lambda_.shape[2-1]
    
    h = lambda_.shape[3-1] - 1
    
    gamma = np.zeros((n,h))
    for i in np.arange(1,n+1).reshape(-1):
        for j in np.arange(1,h+1).reshape(-1):
            # compute gamma[i,j]
            Rvi = 0 * R[1,i]
            for v in np.arange(1,u+1).reshape(-1):
                Rvi = Rvi + R[v,i]
            Pij = np.array([1 + j])
            pr_j = par(Rvi,1 + j)
            while not len(pr_j)==0 :

                Pij = np.array([pr_j,Pij])
                pr_j = par(Rvi,pr_j)

            if len(Pij)==0:
                gamma[i,j] = 0
            else:
                gamma[i,j] = 1
                for li in np.arange(2,len(Pij)+1).reshape(-1):
                    y = 0
                    l_1 = Pij(li - 1)
                    l = Pij(li)
                    for v in np.arange(1,u+1).reshape(-1):
                        for t in np.arange(1,l_1+1).reshape(-1):
                            y = y + lambda_(v,i,t) * R[v,i](t,l)
                    gamma[i,j] = scipy.special.gamma[i,j] * y
    
    return gamma,u,n,h
        
def par(R = None,j = None): 
    # finds the parent of j according to the access probabilities in R
    parent = find(R(np.arange(1,(j - 1)+1),j))
    if len(parent) > 1:
        line_error(mfilename,'A cache has a list with more than one parent, but the structure must be a tree.')
    
    return parent
    
def cache_miss(gamma = None,m = None,lambda_ = None): 
    # M: global miss rate
# MU: per-user miss rate
# MI: per-item miss rate
# pi0: per-item miss probability
    ma = m
    ma[1] = ma(1) + 1
    M = cache_erec(gamma,ma) / cache_erec(gamma,m)
    if len(varargin) > 2:
        u = lambda_.shape[1-1]
        n = lambda_.shape[2-1]
        MU = np.zeros((u,1))
        for v in np.arange(1,u+1).reshape(-1):
            for k in np.arange(1,n+1).reshape(-1):
                pi0[k] = cache_erec(scipy.special.gamma(setdiff(np.arange(1,n+1),k),:),m) / cache_erec(gamma,m)
                MU[v] = MU(v) + (lambda_(v,k,1)) * pi0[k]
        MI = np.zeros((n,1))
        for k in np.arange(1,n+1).reshape(-1):
            MI[k] = MI[k] + sum(lambda_(:,k,1)) * pi0[k]
    
    return M,MU,MI,pi0
    
def cache_miss_asy(gamma = None,m = None,lambda_ = None): 
    # FPI method
    n,h = gamma.shape
    xi = cache_xi_fp(gamma,m)
    MI = np.zeros((n,1))
    for i in np.arange(1,n+1).reshape(-1):
        MI[i] = sum(lambda_(:,i,1)) / (1 + scipy.special.gamma[i,:] * xi)
    
    M = sum(MI)
    if len(varargin) > 2:
        u = lambda_.shape[1-1]
        n = lambda_.shape[2-1]
        MU = np.zeros((u,1))
        for i in np.arange(1,n+1).reshape(-1):
            pi0[i] = 1 / (1 + scipy.special.gamma[i,:] * xi)
            for v in np.arange(1,u+1).reshape(-1):
                MU[v] = MU(v) + lambda_(v,i,1) * pi0[i]
    
    return M,MU,MI,pi0
    
def cache_miss_rayint(gamma = None,m = None,lambda_ = None): 
    # M: global miss rate
# MU: per-user miss rate
# MI: per-item miss rate
    ma = m
    ma[1] = ma(1) + 1
    __,lE,xi = cache_rayint(gamma,m)
    __,lEa,__ = cache_rayint(gamma,ma)
    M = np.exp(lEa - lE)
    if len(varargin) > 2:
        u = lambda_.shape[1-1]
        n = lambda_.shape[2-1]
        pi0 = np.zeros((1,n))
        # compute MU
        MU = np.zeros((u,1))
        if nargout > 1:
            for k in np.arange(1,n+1).reshape(-1):
                if sum(scipy.special.gamma[k,:]) > 0:
                    __,lE1[k] = cache_rayint(scipy.special.gamma(setdiff(np.arange(1,n+1),k),:),m,xi)
                    pi0[k] = np.exp(lE1[k] - lE)
                    if pi0[k] > 1 or pi0[k] < 0:
                        __,lE1[k] = cache_rayint(scipy.special.gamma(setdiff(np.arange(1,n+1),k),:),m)
                        pi0[k] = np.exp(lE1[k] - lE)
                    for v in np.arange(1,u+1).reshape(-1):
                        MU[v] = MU(v) + (lambda_(v,k,1)) * pi0[k]
        # compute MI
        if nargout > 2:
            MI = np.zeros((n,1))
            for k in np.arange(1,n+1).reshape(-1):
                if sum(scipy.special.gamma[k,:]) > 0:
                    MI[k] = MI[k] + sum(lambda_(:,k,1)) * np.exp(lE1[k] - lE)
                else:
                    MI[k] = 0
    
    return M,MU,MI,pi0,lE
    
def cache_mva(gamma = None,m = None): 
    n,h = gamma.shape
    SS = []
    for l in np.arange(1,h+1).reshape(-1):
        SS = ssg_decorate(SS,np.transpose(np.array([np.arange(1,(m(l) + 1)+1)])))
    
    SS = SS - 1
    np.pi = np.zeros((SS.shape[1-1],n))
    pij = np.zeros((SS.shape[1-1],n,h))
    x = np.zeros((1,h))
    E = 1
    #Ecur=SS(1,:);
    for s in np.arange(1,SS.shape[1-1]+1).reshape(-1):
        mcur = SS(s,:)
        for l in np.arange(1,h+1).reshape(-1):
            mcur_l = oner(mcur,l)
            s_l = matchrow(SS,mcur_l)
            if s_l > 0:
                x[l] = mcur(l) / (np.transpose(scipy.special.gamma(:,l)) * np.transpose((1 - np.pi(s_l,:))))
                pij[s,:,l] = np.multiply(np.transpose(scipy.special.gamma(:,l)),(1 - np.pi(s_l,:))) * x(l)
                np.pi[s,:] = np.pi(s,:) + pij(s,:,l)
    
    s = matchrow(SS,m)
    np.pi = np.transpose(np.pi(s,:))
    pij = reshape(pij(s,:,:),n,h)
    pi0 = 1 - np.pi
    if nargout > 2:
        for l in np.arange(1,h+1).reshape(-1):
            for k in np.arange(1,n+1).reshape(-1):
                u[k,l] = x(l) * scipy.special.gamma(k,l)
    
    return np.pi,pi0,pij,x,u,E
    
    
def ssg_decorate(SS = None,SS2 = None): 
    # SS = ssg_decorate(SS1, SS2)
# INPUT:
# SS1 : a state space (n1,...,nk)
# SS2 : a state space (m1,...,mk)
# OUTPUT:
# SS  : a state space (n1,...,nk,m1,...,mk)
# EXAMPLE:
# SS = ssg_closed_single(3,2)
# RR=ssg_renv(10)
# ssg_decorate(SS,RR)
    
    # Copyright (c) 2012-2014, Imperial College London
# All rights reserved.
    
    if len(SS)==0:
        SS = SS2
        return SS
    else:
        if len(SS2)==0:
            return SS
        n1 = SS.shape[1-1]
        m1 = SS.shape[2-1]
        n2 = SS2.shape[1-1]
        m2 = SS2.shape[2-1]
        SS = np.matlib.repmat(SS,n2,1)
        curStates = np.arange(1,n1+1)
        for s in np.arange(1,n2+1).reshape(-1):
            SS[curStates,np.arange[[m1 + 1],[m1 + m2]+1]] = np.matlib.repmat(SS2(s,:),len(curStates),1)
            curStates = (curStates) + n1
    return SS
    
def cache_mva_miss(p = None,m = None,R = None): 
    n = len(p)
    h = len(m)
    if sum(m) == 0 or np.amin(m) < 0:
        Mk = np.ones((1,n))
        M = p * np.transpose(Mk)
        return M,Mk
    
    for j in np.arange(1,h+1).reshape(-1):
        __,Mj = cache_mva_miss(p,oner(m,j),R)
        for k in np.arange(1,n+1).reshape(-1):
            w[k,j] = np.prod(R(np.arange(1,j+1),k)) * p[k] ** j * np.abs(Mj[k])
    
    for j in np.arange(1,h+1).reshape(-1):
        x[j] = 1 / sum(np.abs(w(:,j)))
    
    for k in np.arange(1,n+1).reshape(-1):
        Mk[k] = 1
        for j in np.arange(1,h+1).reshape(-1):
            Mk[k] = Mk[k] - x[j] * m[j] * w[k,j]
    
    Mk = np.abs(Mk)
    M = p * np.transpose(Mk)
    return M,Mk
    
def cache_prob_asy(gamma = None,m = None): 
    # FPI method
    n,h = gamma.shape
    xi = cache_xi_fp(gamma,m)
    prob = np.zeros((n,h))
    for i in np.arange(1,n+1).reshape(-1):
        prob[i,1] = 1 / (1 + scipy.special.gamma[i,:] * xi)
        prob[i,np.arange[2,[1 + h]+1]] = scipy.special.gamma[i,:] * xi / (1 + scipy.special.gamma[i,:] * xi)
    
    return prob
    
def cache_prob_erec(gamma = None,m = None): 
    n,h = gamma.shape
    E = cache_erec(gamma,m)
    for i in np.arange(1,n+1).reshape(-1):
        for j in np.arange(1,h+1).reshape(-1):
            Ei = cache_erec(scipy.special.gamma(setdiff(np.arange(1,n+1),i),:),oner(m,j))
            prob[i,1 + j] = m[j] * scipy.special.gamma[i,j] * Ei / E
        prob[i,1] = np.abs(1 - sum(prob(i,np.arange(2,-1))))
    
    return prob
    
def cache_prob_rayint(gamma = None,m = None,lE = None): 
    n,h = gamma.shape
    if len(varargin) < 3:
        __,lE = cache_rayint(gamma,m)
    
    for i in np.arange(1,n+1).reshape(-1):
        for j in np.arange(1,h+1).reshape(-1):
            __,lEi = cache_rayint(scipy.special.gamma(setdiff(np.arange(1,n+1),i),:),oner(m,j))
            prob[i,1 + j] = m[j] * scipy.special.gamma[i,j] * np.exp(lEi - lE)
        prob[i,1] = np.abs(1 - sum(prob(i,np.arange(2,-1))))
    
    return prob
    
def cache_rayint(gamma = None,m = None,xi0 = None): 
    gamma = scipy.special.gamma(find(np.sum(gamma, 2-1) > 0),:)
    h = len(m)
    n = len(gamma)
    mt = sum(m)
    if n == mt:
        line_warning(mfilename,'The number of items equals the cache capacity')
    
    if len(varargin) < 3:
        xi = cache_xi_bvh(gamma,m)
    else:
        xi = cache_xi_bvh(gamma,m,xi0)
    
    for k in np.arange(1,n+1).reshape(-1):
        S[k] = 0
        for l in np.arange(1,h+1).reshape(-1):
            S[k] = S[k] + scipy.special.gamma(k,l) * xi(l)
    
    ## phi
    phi = 0
    for k in np.arange(1,n+1).reshape(-1):
        phi = phi + np.log(1 + S[k])
    
    phi = phi - np.log(xi) * np.transpose(m)
    ## A
    delta = np.eye(h)
    for j in np.arange(1,h+1).reshape(-1):
        for l in np.arange(1,h+1).reshape(-1):
            C1 = 0
            for k in np.arange(1,n+1).reshape(-1):
                C1 = C1 + scipy.special.gamma[k,j] / (1 + S[k])
            C2 = 0
            for k in np.arange(1,n+1).reshape(-1):
                C2 = C2 + scipy.special.gamma[k,j] * scipy.special.gamma(k,l) / (1 + S[k]) ** 2
            C[j,l] = delta(j,l) * C1 - xi[j] * C2
    
    ##
    Z = np.exp(phi) * np.sqrt(2 * np.pi) ** (- h) * np.prod(factorial(m)) / np.prod(np.sqrt(xi)) / np.sqrt(det(C))
    lZ = (- h) * np.log(np.sqrt(2 * np.pi)) + (phi) + sum(factln((m))) - sum(np.log(np.sqrt(xi))) - np.log(np.sqrt(det(C)))
    lZ = real(lZ)
    
    if not isfinite(lZ) :
        #keyboard
        pass
    
    return Z,lZ,xi
    
def cache_xi_bvh(gamma = None,m = None,tmax = None): 
    if len(varargin) < 3:
        tmax = Inf
    
    n = gamma.shape[1-1]
    f = m / n
    h = f.shape[2-1]
    pp = np.zeros((h + 1,n))
    pp[1,:] = np.ones((1,n))
    for i in np.arange(1,h+1).reshape(-1):
        pp[i + 1,:] = scipy.special.gamma(:,i)
    
    z_old = np.zeros((1,h + 1))
    z = np.ones((1,h + 1))
    T = tic
    while (np.amax(np.abs(z - z_old) > 10 ** (- 12) * np.amax(np.abs(z_old)))):

        if toc(T) > tmax:
            break
        z_old = z
        temp = n * z * pp
        for i in np.arange(1,h+1).reshape(-1):
            a = temp - n * z(i + 1) * pp(i + 1,:)
            Fi = sum(pp(i + 1,:) / (n * pp(i + 1,:) + a))
            if (Fi > f[i]):
                zi_min = 0
                zi_max = 1
            else:
                zi_min = 1
                zi_max = 2
                while (sum(zi_max * pp(i + 1,:) / (n * zi_max * pp(i + 1,:) + a)) < f[i]):

                    zi_min = zi_max
                    zi_max = zi_max * 2

            for x in np.arange(1,50+1).reshape(-1):
                z[i + 1] = (zi_min + zi_max) / 2
                if (sum(z(i + 1) * pp(i + 1,:) / (n * z(i + 1) * pp(i + 1,:) + a)) < f[i]):
                    zi_min = z(i + 1)
                else:
                    zi_max = z(i + 1)

    
    z = z(np.arange(2,-1))
    return z
        
def cache_xi_fp(gamma = None,m = None,xi = None): 
    n,h = gamma.shape
    tol = 1e-14
    pi0 = np.ones((1,n)) / (h + 1)
    pij = np.zeros((n,h))
    xi = np.zeros((1,h))
    if len(varargin) < 3:
        for l in np.arange(1,h+1).reshape(-1):
            xi[l] = m(l) / mean(scipy.special.gamma(:,l)) / (n + sum(m) - 1)
    
    for it in np.arange(1,10000.0+1).reshape(-1):
        pi0_1 = pi0
        xi = m / (pi0_1 * gamma)
        pij = np.abs(np.multiply(gamma,np.matlib.repmat(xi,n,1))) / np.abs(1 + gamma * np.transpose(xi))
        pi0 = np.amax(tol,1 - np.transpose(np.sum(pij, 2-1)))
        DELTA = norm(np.abs(1 - pi0 / pi0_1),1)
        if DELTA < tol:
            xi[xi < 0] = tol
            return xi,pi0,pij,it
    
    xi[xi < 0] = tol
    return xi,pi0,pij,it
    