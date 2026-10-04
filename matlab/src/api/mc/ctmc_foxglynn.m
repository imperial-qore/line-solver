function [pi,left,right,w]=ctmc_foxglynn(pi0,Q,t,tol,maxiter)
% [PI,LEFT,RIGHT,W]=CTMC_FOXGLYNN(PI0,Q,T,TOL,MAXITER)
%
% Transient distribution of the CTMC by uniformization, with Poisson weights
% and truncation points obtained by the Fox-Glynn algorithm, with Jansen's
% correction to the right tail estimate. Neither exp(-q*t) nor (q*t)^k/k! is
% ever formed, so the method is free of overflow and underflow and needs no
% horizon splitting. MAXITER caps the right truncation point; pass a
% nonpositive value (or omit it) to leave it uncapped.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin<4%~exist('tol','var')
    tol = 1e-12;
end
if nargin<5%~exist('maxiter','var')
    maxiter = -1;
end
if tol<=0
    tol = 1e-12;
end
q=1.1*max(abs(diag(Q)));
lambda=q*t;
if q<=0 || lambda<=0
    pi = pi0;
    left = 0;
    right = 0;
    w = 1;
    return
end
[left,right,w] = foxglynn_weights(lambda,tol,maxiter);
Qs=speye(size(Q))+sparse(Q)/q;
pi=zeros(size(pi0));
P=pi0;
for k=0:right
    if k>=left
        pi = pi + w(k-left+1)*P;
    end
    if k<right
        P = P*Qs;
    end
end
end
