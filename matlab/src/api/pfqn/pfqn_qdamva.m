%{
%{
 % @file pfqn_qdamva.m
 % @brief QD-AMVA: queue-dependent approximate mean value analysis.
%}
%}

%{
%{
 % @brief QD-AMVA, the queue-dependent AMVA of Casale-Perez-Wang (IFIP
 %        PERFORMANCE 2015), on a closed multiclass product-form network.
 % @details
 % Schweitzer/Bard core in which the class-r demand at station k is scaled by
 % the queue-dependence term g_k evaluated at the arrival-instant total queue
 % length, g = pfqn_lldfun(1 + delta*sum(Q(k,:)), mu).
 %
 % SETTING MU TO A CONSTANT ROW RECOVERS PLAIN SCHWEITZER AMVA ONLY FOR A SINGLE
 % CLASS. pfqn_lldfun does skip a constant row, so g == 1 there, but the residence
 % time that remains is 1 + delta*sum(Q(k,:)) with ONE aggregate delta =
 % (sum(N)-1)/sum(N) applied to the whole arrival-instant queue, where Bard-
 % Schweitzer shrinks the TAGGED class alone:
 %   1 + sum_{s~=r} Q(k,s) + (N(r)-1)/N(r) * Q(k,r).
 % The two coincide iff K == 1. Measured over 40 random three-class instances,
 % pfqn_qdamva(L,N,Z,ones) departs from pfqn_bs by up to 0.217 in absolute queue
 % length, and is the LESS accurate of the two on single-server multiclass models
 % (mean relative error on Q 0.069 against 0.056 at R = 3), the aggregate delta
 % buying nothing once g == 1. This is the QD-AMVA closure, not a defect of the
 % implementation, but do not use the function as a Schweitzer oracle for K > 1.
 %
 % MU IS A DIMENSIONLESS RATE MULTIPLIER, NOT A RATE. mu(k,n) is the factor by
 % which station k serves faster when it holds n jobs. Two traps follow from
 % pfqn_lldfun:
 %   - it SKIPS a station whose mu row is constant (its `range(...)>0` gate), so
 %     a single-server station must be ones(1,smax) and a c-server station
 %     min(1:smax, c). Passing a c-server station a constant row silently
 %     returns g=1, i.e. a single server.
 %   - smax = size(mu,2) must be at least ceil(sum(N)) or its interp1 clamps the
 %     population and the top of the rate curve is never reached.
 %
 % Delay stations are carried in Z, not as rows of L. Closed classes only: an
 % infinite N(r) is not supported.
 %
 % @fn pfqn_qdamva(L, N, Z, mu, Q0, tol, maxiter)
 % @param L (M x R) service demand matrix.
 % @param N (1 x R) population vector, finite.
 % @param Z (1 x R) think time vector.
 % @param mu (M x smax) queue-dependent rate multipliers.
 % @param Q0 (M x R) initial guess for the queue lengths.
 % @param tol Convergence tolerance on the queue lengths (default 1e-6).
 % @param maxiter Maximum number of iterations (default 1e4).
 % @return Q (M x R) mean queue lengths.
 % @return X (1 x R) per-class throughputs.
 % @return U (M x R) per-class utilizations, carrying the g scaling.
 % @return iter Number of iterations performed.
 % @return R (M x R) per-class residence times, Q = X.*R.
%}
%}
function [Q,X,U,iter,R] = pfqn_qdamva(L,N,Z,mu,Q0,tol,maxiter)
% [Q,X,U,ITER,R] = PFQN_QDAMVA(L,N,Z,MU,Q0,TOL,MAXITER)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,K] = size(L);
N = N(:)';
if nargin < 3 || isempty(Z)
    Z = zeros(1,K);
end
Z = Z(:)';
if nargin < 4
    mu = [];
end
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
if nargin < 7 || isempty(maxiter)
    maxiter = 1e4;
end

if nargin < 5 || isempty(Q0)
    % Ltot=0 for a class with no demand anywhere: L./Ltot is NaN and the
    % iteration never recovers. Such a column arises routinely in a layered
    % fixed point, where a caller can start with no work at the layer station.
    Ltot = sum(L,1);
    Q = zeros(M,K);
    nz = Ltot > 0;
    if any(nz)
        Q(:,nz) = L(:,nz) ./ repmat(Ltot(nz),M,1) .* repmat(N(nz),M,1);
    end
else
    Q = Q0;
end

X = zeros(1,K);
U = zeros(M,K);
R = zeros(M,K);
iter = 0;
if sum(N) <= 0
    Q = zeros(M,K); % delta is undefined on an empty population
    return
end

delta = (sum(N) - 1) / sum(N);

% Q*10 as the sentinel, as in the reference, stalls on an all-zero seed: the
% loop would exit before its first pass. Offset instead.
Q_1 = Q + 10*(1+tol);
while max(max(abs(Q-Q_1))) > tol && iter < maxiter
    iter = iter + 1;
    Q_1 = Q;

    % arrival-instant total queue length, class-independent
    Ak = 1 + delta * sum(Q,2);
    g = pfqn_lldfun(Ak, mu);

    for r = 1:K
        for k = 1:M
            R(k,r) = L(k,r) * g(k) * (1 + delta * sum(Q(k,:)));
        end
        denom = Z(r) + sum(R(:,r));
        if denom > 0
            X(r) = N(r) / denom;
        else
            X(r) = 0;
        end
        for k = 1:M
            Q(k,r) = X(r) * R(k,r);
            U(k,r) = L(k,r) * g(k) * X(r);
        end
    end
end
end
