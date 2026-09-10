%{
%{
 % @file pfqn_busyp.m
 % @brief Mean busy period of order n for a subnetwork of a product-form network.
%}
%}

%{
%{
 % @brief Mean busy period of order n for a subnetwork of a product-form network.
 % @fn pfqn_busyp(alpha, mu, P, N, subnet, n, gamma, options)
 % @param alpha Relative arrival rates (1xJ), solution of the traffic equations.
 % @param mu Load-dependent service rates (JxK), mu(j,k) with k jobs at node j,
 %           or a handle mu(j,kvec) when the rates do not saturate.
 % @param P Routing matrix (JxJ).
 % @param N Population (scalar, Inf for an open network).
 % @param subnet Indexes of the nodes forming the subnetwork.
 % @param n Busy period order(s), 1 <= n <= N.
 % @param gamma External arrival rates (1xJ), empty for a closed network.
 % @param tol Relative tolerance of the open-network tail truncation.
 % @return b Mean busy period duration(s), same size as n.
 % @return lG Log normalizing constants of the subnetwork, orders 0..K.
 % @return lH Log normalizing constants of the complement, orders 0..N.
%}
%}
function [b, lG, lH] = pfqn_busyp(alpha, mu, P, N, subnet, n, gamma, tol)
% [B,LG,LH] = PFQN_BUSYP(ALPHA, MU, P, N, SUBNET, N, GAMMA, TOL)
%
% Mean duration of the busy period of order n for the subnetwork SUBNET, that
% is the time from the instant a job entering the subnetwork finds n-1 jobs in
% it up to the next instant when fewer than n jobs remain in it.
%
% Implements H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks:
% Mean Value Analysis", J. ACM 35(3), 1988, Theorem 1 (closed Gordon-Newell
% network) and Theorem 3 (open Jackson network). Both are evaluated in the log
% domain, which serves the same purpose as the ratio recursions of Corollaries
% 2 and 4, namely avoiding the overflow of the individual normalizing constants.
%
% The paper is single-chain: ALPHA is the stochastic solution of x*P = x for a
% closed network and the solution of x = GAMMA + x*P for an open one, and every
% node is a state-dependent single-server FCFS station. By the insensitivity
% noted in the paper (Section 5) the result depends on the service processes
% only through the rates MU.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 7
    gamma = [];
end
if nargin < 8 || isempty(tol)
    % NOT options.tol: that is the iterative-solver tolerance (1e-4 by default),
    % which truncates the open-network tail sum three orders too early. The
    % same default is hardcoded in the Python and JAR twins so the three agree.
    tol = 1e-12;
end

alpha = alpha(:)';
J = length(alpha);
if isempty(N)
    N = Inf;
end
isClosed = isfinite(N);

subnet = unique(subnet(:)');
if isempty(subnet)
    line_error(mfilename, 'The subnetwork must be non-empty.');
end
if isClosed && numel(subnet) >= J
    % a closed network needs jobs outside the subnetwork to start a busy period
    line_error(mfilename, 'In a closed network the subnetwork must be a proper subset of the nodes.');
end
if any(subnet < 1) || any(subnet > J)
    line_error(mfilename, 'The subnetwork indexes are out of range.');
end
compl = setdiff(1:J, subnet);

if isClosed
    if N ~= round(N) || N < 1
        line_error(mfilename, 'The population of a closed network must be a positive integer.');
    end
    if any(n < 1) || any(n > N) || any(n ~= round(n))
        line_error(mfilename, 'The busy period order must be an integer in 1..N.');
    end
else
    if isempty(gamma)
        line_error(mfilename, 'An open network requires the external arrival rates gamma.');
    end
    if any(n < 1) || any(n ~= round(n))
        line_error(mfilename, 'The busy period order must be a positive integer.');
    end
end

% A(I) for a closed network, C(I) for an open one: both are the total rate at
% which jobs enter the subnetwork from outside it, which is what starts a busy
% period. The closed network has no external stream, so gamma contributes 0.
inflow = sum(sum(alpha(compl)' .* P(compl, subnet)));
if ~isempty(gamma)
    gamma = gamma(:)';
    inflow = inflow + sum(gamma(subnet));
end
if inflow <= 0
    line_error(mfilename, 'No job ever enters the subnetwork, its busy period is undefined.');
end

if isClosed
    lG = local_lgvec(alpha(subnet), local_rates(mu, subnet, N), N);
    lH = local_lgvec(alpha(compl), local_rates(mu, compl, N), N);
    b = zeros(size(n));
    for t = 1:numel(n)
        nt = n(t);
        % Theorem 1: sum_{m=n}^{N} G(m,I) H(N-m,I) over G(n-1,I) H(N-n,I) A(I)
        num = local_lse(lG(n(t) +1 : N +1) + lH(N-nt +1 : -1 : 0 +1));
        b(t) = exp(num - lG(nt-1 +1) - lH(N-nt +1) - log(inflow));
    end
else
    K = local_trunc(alpha(subnet), mu, subnet, max(n), tol);
    lG = local_lgvec(alpha(subnet), local_rates(mu, subnet, K), K);
    lH = [];
    b = zeros(size(n));
    for t = 1:numel(n)
        % Theorem 3: sum_{m=n}^{Inf} G(m,I) over G(n-1,I) C(I), the tail being
        % summed up to the truncation order fixed by local_trunc.
        num = local_lse(lG(n(t) +1 : K +1));
        b(t) = exp(num - lG(n(t)-1 +1) - log(inflow));
    end
end
end

function murows = local_rates(mu, idx, K)
% Rates of the selected nodes for populations 1..K. A rate table shorter than K
% keeps its last rate, the saturated-server convention; a node whose rate does
% not saturate (an infinite server) must be supplied as a handle instead.
if isa(mu, 'function_handle')
    murows = zeros(numel(idx), K);
    for t = 1:numel(idx)
        murows(t, :) = mu(idx(t), 1:K);
    end
    return
end
murows = mu(idx, :);
if size(murows, 2) < K
    murows = [murows, repmat(murows(:, end), 1, K - size(murows, 2))];
end
murows = murows(:, 1:K);
end

function lg = local_lgvec(alpha, mu, K)
% lg(m+1) = log sum over the compositions n_1+...+n_L = m of the product over
% the nodes of prod_{k=1}^{n_i} alpha_i/mu_i(k), the G(m,I) and H(m,I) of the
% paper. The nodes are convolved one at a time in the log domain.
L = length(alpha);
lg = [0, -Inf(1, K)];
for i = 1:L
    % per-node sequence log prod_{k=1}^{m} alpha_i/mu_i(k)
    li = [0, cumsum(log(alpha(i)) - log(mu(i, 1:K)))];
    lgnew = -Inf(1, K+1);
    for m = 0:K
        lgnew(m +1) = local_lse(lg(m +1 : -1 : 0 +1) + li(0 +1 : m +1));
    end
    lg = lgnew;
end
end

function K = local_trunc(alpha, mu, subnet, nmax, tol)
% Truncation order of the open-network sum sum_{m>=n} G(m,I). The subnetwork
% terms decay geometrically at the rate of its most utilized node once the
% rates saturate, so the truncation is grown until the geometric tail estimate
% is negligible against the partial sum.
K = max(nmax + 8, 16);
murows = local_rates(mu, subnet, K);
rho = alpha ./ murows(:, end)';
if max(rho) >= 1
    line_error(mfilename, 'The subnetwork is not stable, its busy period is infinite.');
end
while true
    lg = local_lgvec(alpha, local_rates(mu, subnet, K), K);
    % decay rate read off the last two orders, the exact ratio for a saturated
    % single-server subnetwork and an upper estimate otherwise
    r = exp(lg(K +1) - lg(K-1 +1));
    if ~(r < 1)
        r = max(rho);
    end
    ltail = lg(K +1) + log(r) - log(1 - r);
    if ltail - local_lse(lg(nmax +1 : K +1)) < log(tol)
        break
    end
    K = 2*K;
    if K > 1e6
        line_error(mfilename, 'The open busy period sum did not converge, the subnetwork is nearly saturated.');
    end
end
end

function s = local_lse(v)
% log-sum-exp, stable when every entry is -Inf
v = v(:)';
m = max(v);
if ~isfinite(m)
    s = m;
    return
end
s = m + log(sum(exp(v - m)));
end
