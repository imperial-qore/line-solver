%{
%{
 % @file pfqn_ncdtld.m
 % @brief Normalizing constants of a discrete-time closed cycle of state
 %        dependent Bernoulli servers.
%}
%}

function [lG, G, W, Gc, Wa] = pfqn_ncdtld(P, N)
%{
%{
 % @brief Convolution over the discrete-time product form with load dependence.
 % @fn pfqn_ncdtld(P, N)
 % @param P J-by-N matrix of service probabilities, P(j,n)=p_j(n) in (0,1).
 % @param N Number of customers cycling in the J nodes.
 % @return lG Logarithm of the normalizing constant G(N,J).
 % @return G Row vector of normalizing constants G(k), k=0..N, in a common scale.
 % @return W J-by-(N+1) matrix of time-stationary node weights, same scale.
 % @return Gc J-by-(N+1) matrix of complement constants, cycle without node j.
 % @return Wa J-by-(N+1) matrix of arrival node weights, same scale as W.
%}
%}
% PFQN_NCDTLD evaluates the discrete-time product form of a closed cycle of J
% state dependent Bernoulli servers. With q_j(n)=1-p_j(n) the queue length
% vector has the stationary law of Daduna (2001), theorem 3.2,
%
%   pi(n_1,...,n_J) = prod_j w_j(n_j) / G(N,J),
%   w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h),   w_j(0) = 1,
%
% which reduces to the state independent form of PFQN_NCDT when p_j(n) does
% not depend on n. Note that w_j misses the factor q_j(n_j) in the numerator,
% so the extra 1/q_j on a busy node in the state independent case is not tied
% to the node being non-empty but to its actual queue length.
%
% The constants are assembled by truncated convolution of the per-node weight
% vectors, together with the J complement constants (the cycle with node j
% removed) obtained from a prefix/suffix pass. All of G, W, Gc and Wa are
% returned in one common scale, so that the marginal law
%
%   P(X_j = n) = W(j,n+1) * Gc(j,N-n+1) / G(N+1)
%
% is exact whatever the scale; lG carries the true log constant. Deconvolution
% is never used, so a node with a near-unit service probability does not
% destroy the accuracy of the other marginals.
%
% Wa holds the weights of the arrival distribution of proposition 3.5,
% v_j(n) = w_j(n) q_j(n), the law seen by a customer joining node j with
% himself not counted:
%
%   Parr(X_j = n) = Wa(j,n+1) * Gc(j,N-n) / sum_m Wa(j,m+1) Gc(j,N-m),
%
% defined on n = 0..N-1. It is not the time-stationary law, which is the
% discrete-time departure from the arrival theorem of continuous time.
%
% Examples:
%   P = repmat([0.5;0.25;0.7], 1, 5);
%   [lG, G, W, Gc] = pfqn_ncdtld(P, 5);
%   W(1,2)*Gc(1,5)/G(6)              % P(X_1 = 1)
%
% Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
% Springer 2001, theorem 3.2 and proposition 3.5.
%
% See also PFQN_NCDT, QSYS_BERNOULLI1, PFQN_NCLD
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isnumeric(N) || ~isscalar(N) || ~isreal(N) || N < 0 || floor(N) ~= N
    line_error(mfilename, 'N must be a non-negative integer');
end
if ~isnumeric(P) || ~isreal(P) || isempty(P)
    line_error(mfilename, 'P must be a non-empty real matrix of service probabilities');
end
[J, cols] = size(P);
if cols < N
    line_error(mfilename, 'P must supply p_j(n) for every n = 1..%d, but has %d columns', N, cols);
end
if N == 0
    lG = 0; G = 1; W = ones(J,1); Gc = ones(J,1); Wa = ones(J,1);
    return
end
Pn = P(:, 1:N);
if any(~isfinite(Pn(:))) || any(Pn(:) <= 0) || any(Pn(:) > 1)
    line_error(mfilename, 'service probabilities must be real and in the interval (0,1]');
end
if N >= 1 && any(any(Pn(:, 1:N-1) >= 1))
    % p_j(n)=1 is admissible only at the last reachable population, where the
    % missing q_j(n) never multiplies any weight.
    line_error(mfilename, 'service probabilities below the population bound must be strictly less than 1');
end

% Per-node log weights of theorem 3.2, then a per-node shift so that every
% weight vector peaks at one. The shifts cancel in every ratio below.
lw = zeros(J, N+1);
for j = 1:J
    cq = 0; cp = 0;
    for n = 1:N
        cp = cp + log(Pn(j, n));
        lw(j, n+1) = cq - cp;
        cq = cq + log(1 - Pn(j, n));
    end
end
shift = max(lw, [], 2);
W = exp(lw - shift(:, ones(1, N+1)));

% Arrival weights v_j(n) = w_j(n) q_j(n) of proposition 3.5.
Wa = W;
for j = 1:J
    for n = 1:N
        Wa(j, n+1) = W(j, n+1) * (1 - Pn(j, n));
    end
end

% Prefix/suffix convolution: pre(j,:) covers nodes 1..j-1, suf(j,:) covers
% nodes j+1..J, so the complement of node j is their convolution and the full
% constant is the convolution of pre(J+1,:) with node J.
pre = zeros(J+1, N+1); pre(1, 1) = 1;
for j = 1:J
    pre(j+1, :) = dt_conv(pre(j, :), W(j, :), N);
end
suf = zeros(J+1, N+1); suf(J+1, 1) = 1;
for j = J:-1:1
    suf(j, :) = dt_conv(suf(j+1, :), W(j, :), N);
end
G = pre(J+1, :);
Gc = zeros(J, N+1);
for j = 1:J
    Gc(j, :) = dt_conv(pre(j, :), suf(j+1, :), N);
end

lG = log(G(N+1)) + sum(shift);
end

% ==========================================================================
function c = dt_conv(a, b, N)
% Convolution of two population tables truncated at N.
c = zeros(1, N+1);
for k = 0:N
    s = 0;
    for m = 0:k
        s = s + a(m+1) * b(k-m+1);
    end
    c(k+1) = s;
end
end
