%{
%{
 % @file pfqn_bklc.m
 % @brief Birman-Kogan load concealment algorithm: multichain solved as single chain problems.
%}
%}

%{
%{
 % @brief Birman-Kogan load concealment algorithm: multichain solved as single chain problems.
 %
 % Birman and Kogan (Stochastic Models 8(3):543-563, 1992), Algorithm 2. The
 % saddle point analysis of Corollary 2 shows that chain l may be solved on
 % its own provided every station is slowed by the residual capacity the other
 % chains leave it, A_i = 1 - sum_{k != l} L(i,k)*X_k. Sweeping the chains in
 % Gauss-Seidel order and iterating to a fixed point is the load concealment
 % algorithm; the paper's contribution is the asymptotic argument
 % that says when it is exact, and the extension to state dependent servers.
 %
 % The single chain subproblem is solved either exactly (MVA) or by the
 % uniform expansion of PFQN_BKUE.
 %
 % @fn pfqn_bklc(L, N, Z, method, tol, maxiter)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector (1 x classes).
 % @param Z Think time vector (default: zeros).
 % @param method Single chain solver, 'mva' (default) or 'ue'.
 % @param tol Convergence tolerance on the throughputs (default: 1e-10).
 % @param maxiter Maximum number of sweeps (default: 1000).
 % @return X Chain throughputs.
 % @return Q Mean queue lengths (stations x classes).
 % @return U Utilizations (stations x classes).
 % @return it Number of sweeps performed.
%}
%}
function [X,Q,U,it] = pfqn_bklc(L,N,Z,method,tol,maxiter)
% [X,Q,U,IT] = PFQN_BKLC(L,N,Z,METHOD,TOL,MAXITER)
%
% Load concealment reduction of a multichain closed network.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,R] = size(L);
if nargin < 3 || isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if nargin < 4 || isempty(method)
    method = 'mva';
end
if nargin < 5 || isempty(tol)
    tol = 1e-10;
end
if nargin < 6 || isempty(maxiter)
    maxiter = 1000;
end
X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R); it = 0;
if isempty(L) || sum(N) == 0
    return
end
% Step 1: the saddle point utilizations of Corollary 1 seed the iteration
try
    [~,~,X] = pfqn_bk(L,N,Z);
catch
    X = zeros(1,R);
end
ok = isfinite(X) & X >= 0;
X(~ok) = 0;
for r = 1:R
    if X(r) == 0 && N(r) > 0
        X(r) = N(r)/(Z(r) + sum(L(:,r)));
    end
end
% A chain cannot draw more than the capacity of its own slowest station
for r = 1:R
    cap = max(L(:,r));
    if cap > 0
        X(r) = min(X(r), 1/cap);
    end
end

for it = 1:maxiter
    Xold = X;
    for l = 1:R
        if N(l) == 0
            X(l) = 0; Q(:,l) = 0;
            continue
        end
        % Step 2a: residual capacity left to chain l at every station (eq. 32-33)
        A = 1 - (L*X' - L(:,l)*X(l));
        A = max(A, GlobalConstants.FineTol);
        D = L(:,l)./A;
        % Step 2b: solve the single chain network with the concealed rates
        switch method
            case 'ue'
                [Xl,Ql] = local_ue(D,N(l),Z(l));
            otherwise
                [Xl,Ql] = pfqn_mva(D,N(l),Z(l));
        end
        X(l) = Xl;
        Q(:,l) = Ql(:);
    end
    if max(abs(X - Xold)) <= tol*max(1,max(abs(X)))
        break
    end
end
U = L .* repmat(X,M,1);
end

function [X,Q] = local_ue(D,N,Z)
% throughputs from the uniform expansion at successive populations, then the
% queue lengths from the mean value recursion they imply
M = numel(D);
Q = zeros(M,1);
X = 0;
lGprev = 0;
for n = 1:N
    [~,lGn] = pfqn_bkue(D,n,Z);
    X = exp(lGprev - lGn);
    Q = D(:).*X.*(1 + Q);
    lGprev = lGn;
end
end
