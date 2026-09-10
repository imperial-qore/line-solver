%{
%{
 % @file pfqn_scb.m
 % @brief Dowdy-Carlson-Krantz-Tripathi (1992) single-class bounds on the
 %        performance of the multiclass system a single-class model aggregates.
%}
%}

function [Xlo,Xhi,Ulo,Uhi] = pfqn_scb(L,N)
%{
%{
 % @brief Bracket on the total throughput and on the per-device utilizations of
 %        the UNKNOWN multiclass system whose single-class counterpart has
 %        demand vector L at population N. L. W. Dowdy, B. M. Carlson, A. T.
 %        Krantz, S. K. Tripathi, "Single-Class Bounds of Multi-Class Queuing
 %        Networks", J. ACM 39(1):188-213, 1992.
 %
 %        SEMANTICS DIFFER FROM EVERY OTHER pfqn_* BOUND. aba/bjb/gb/... bracket
 %        the exact solution OF THE GIVEN MODEL; this brackets the multiclass
 %        system that the given single-class model aggregates. The lower side is
 %        therefore the EXACT single-class solution, not an approximation of it.
 %
 %        Theorem 2 / Corollary 2: aggregating an R-class model into its
 %        single-class counterpart (demands weighted by the relative class
 %        throughputs) can only understate performance, U_k,1 <= U_k,R and
 %        X_1 <= X_R, and Corollary 1 makes the utilization ratio uniform,
 %        U_k,R/U_k,1 = X_R/X_1 for every k. Theorem 3 (their Expression 3)
 %        caps the relative throughput error at (m-1)/(N+m-1), m = min(N,K),
 %        independently of the demands, hence X_R <= X_1*(N+m-1)/N. The single
 %        server capacity U_k,R <= 1 caps the same ratio at 1/(X_1*max(L)), and
 %        that cap is tight on the paper's own worst case (m saturated devices,
 %        where D_k,1 = 1/X_R for every k), so both are applied.
 % @fn pfqn_scb(L, N)
 % @param L Service demand vector of the single-class model (K x 1), queueing
 %          stations only. Delay stations are not admitted: Theorem 3 rests on
 %          the delay-free balanced-network throughput N/((N+m-1)D).
 % @param N Total population (scalar, N >= 1).
 % @return Xlo Lower bound on multiclass total throughput X_R (= exact X_1).
 % @return Xhi Upper bound on X_R.
 % @return Ulo Lower bound on the per-device utilizations U_k,R (K x 1).
 % @return Uhi Upper bound on U_k,R (K x 1).
%}
%}
L = L(:);
K = numel(L);
if K == 0
    line_error(mfilename,'pfqn_scb requires at least one queueing station.');
end
N = round(N);
if N < 1
    line_error(mfilename,'pfqn_scb requires N >= 1.');
end

% Exact single-class MVA at Z=0. This IS the lower bound (Theorem 2), so it is
% computed exactly rather than bounded: a bounded X1 would not bracket X_R.
Q = zeros(K,1);
X1 = 0;
for n = 1:N
    Rk = L .* (1 + Q);
    X1 = n / sum(Rk);
    Q = X1 * Rk;
end
U1 = X1 * L;

m = min(N,K);
ratio = (N + m - 1) / N;             % Theorem 3, Expression (3)
Dmax = max(L);
if X1 * Dmax > 0
    % U_k,R <= 1 with the uniform ratio of Corollary 1. Tight at the worst case.
    ratio = min(ratio, 1 / (X1 * Dmax));
end

Xlo = X1;
Xhi = X1 * ratio;
Ulo = U1;
Uhi = U1 * ratio;
end
