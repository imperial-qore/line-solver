%{
%{
 % @file pfqn_busyp_clw.m
 % @brief Busy period of a subnetwork from point evaluations of the normalizing
 %        constant (Choudhury-Leung-Whitt inversion).
%}
%}

%{
%{
 % @brief Busy period of a subnetwork through normalizing-constant evaluations.
 % @fn pfqn_busyp_clw(alpha, mu, P, N, subnet, n, gamma, isdelay, method)
 % @param alpha Relative arrival rates (JxR), one column per chain.
 % @param mu Service rates (JxR), the chain-r rate at node j.
 % @param P Routing matrix (JxJ), or a 1xR cell of per-chain matrices.
 % @param N Population per chain (1xR), Inf entries for an open chain.
 % @param subnet Indexes of the nodes forming the subnetwork.
 % @param n Busy period order(s).
 % @param gamma External arrival rates (JxR), empty for a closed network.
 % @param isdelay Logical (Jx1), true at an infinite server.
 % @param method Normalizing-constant method of the point evaluations ('clw').
 % @return b Mean busy period duration(s), same size as n.
%}
%}
function b = pfqn_busyp_clw(alpha, mu, P, N, subnet, n, gamma, isdelay, method)
% B = PFQN_BUSYP_CLW(ALPHA, MU, P, N, SUBNET, N, GAMMA, ISDELAY, METHOD)
%
% WHAT THIS BUYS OVER PFQN_BUSYP / PFQN_BUSYP_MULTICLASS. Those walk the whole
% population ladder (the whole lattice, multichain) because the numerator sums
% over {|m| >= n}. The complement of that set is the SHELLS |m| <= n-1, and
% summing the product form over the WHOLE lattice is the full network's own
% normalizing constant, G_I and H convolving to it:
%
%   sum_{m : |m| >= n} G_I(m) H(N-m) = G(N) - sum_{m : |m| <= n-1} G_I(m) H(N-m)
%
% so the busy period of order n needs only the n lowest shells plus ONE
% evaluation of G(N). The ordinary busy period n=1 collapses to three constants,
%
%   b(1,I) = [G(N) - H(N)] / sum_r A_r(I) H(N-e_r)
%
% all of them point evaluations at or near the full population, which is what
% the normalizing-constant methods are built for. This routine calls CLW
% (Choudhury-Leung-Whitt, J. ACM 42, 1995, numerical inversion of the generating
% function); any method returning lG(N) can take its place. The cost stops
% depending on N: O(shells up to n-1) plus O(n*R) constant evaluations, against
% O(lattice) for the ladder routines.
%
% THE OPEN CASE NEEDS NO INVERSION. The subnetwork's constant sequence has
% generating function g(z) = prod_{i in I} f_i(z), and the tail is g(1) minus a
% partial sum,
%
%   b(n,I) = [g_I(1) - sum_{m=0}^{n-1} G_I(m)] / [G_I(n-1) C(I)],
%
% with f_i(1) = 1/(1-rho_i) at a single server and exp(rho_i) at an infinite
% one. That removes the tail TRUNCATION of the ladder routine, not just its
% cost: the tail is exact.
%
% ACCURACY. The numerator is a difference of two nearly equal quantities when
% the level set is unlikely, so the relative error grows with n: measured 6e-12
% at n=1 against 2.7e-08 at n=N on a three-station closed model at N=20. Cost
% grows with n too, so the routine is most accurate where it is fastest.
%
% SCOPE. CLW's generating function covers single-server and infinite-server
% stations, so a general load-dependent scaling belongs to
% PFQN_BUSYP_MULTICLASS. The identity is for the AGGREGATE level set: a
% per-class one has complement {m_r <= n-1}, the whole lattice in the other
% chains, which buys nothing.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 7
    gamma = [];
end
[J, R] = size(alpha);
if nargin < 8 || isempty(isdelay)
    isdelay = false(J, 1);
end
if nargin < 9 || isempty(method)
    method = 'clw';
end
isdelay = logical(isdelay(:));
N = N(:)';
isClosed = all(isfinite(N));
isOpen = all(~isfinite(N));
if ~isClosed && ~isOpen
    line_error(mfilename, 'A mixed model needs the lattice routine pfqn_busyp_multiclass.');
end

subnet = unique(subnet(:)');
if isempty(subnet)
    line_error(mfilename, 'The subnetwork must be non-empty.');
end
if isClosed && numel(subnet) >= J
    line_error(mfilename, 'In a closed network the subnetwork must be a proper subset of the nodes.');
end
compl = setdiff(1:J, subnet);

% demands L(i,r) = alpha(i,r)/mu(i,r), zero where chain r does not visit node i
L = zeros(J, R);
visited = alpha > 0 & mu > 0;
L(visited) = alpha(visited) ./ mu(visited);

% A_r(I): the chain-r rate at which jobs enter the subnetwork from outside it
A = zeros(1, R);
for r = 1:R
    if iscell(P)
        Pr = P{r};
    else
        Pr = P;
    end
    A(r) = sum(sum(alpha(compl, r) .* Pr(compl, subnet)));
    if ~isempty(gamma)
        A(r) = A(r) + sum(gamma(subnet, r));
    end
end
if sum(A) <= 0
    line_error(mfilename, 'No job ever enters the subnetwork, its busy period is undefined.');
end

if isOpen
    % g_I(1) in closed form, so the tail is exact rather than truncated
    rho = sum(L(subnet, :), 2)';
    queues = ~isdelay(subnet)';
    if any(rho(queues) >= 1)
        line_error(mfilename, 'The subnetwork is not stable, its busy period is infinite.');
    end
    lg1 = 0;
    for t = 1:numel(subnet)
        if isdelay(subnet(t))
            lg1 = lg1 + rho(t);
        else
            lg1 = lg1 - log1p(-rho(t));
        end
    end
    lseq = local_open_coeff(rho, isdelay(subnet)', max(n));
    b = zeros(size(n));
    for t = 1:numel(n)
        head = local_lse(lseq(1:n(t)));
        tail = lg1 + log1p(-exp(head - lg1));
        b(t) = exp(tail - lseq(n(t)) - log(sum(A)));
    end
    return
end

kmax = max(n) - 1;
bound = min(N, max(kmax, 0));
Mvec = local_lattice(bound);
lGlow = local_lgvec(L(subnet, :), isdelay(subnet), Mvec);
level = sum(Mvec, 2)';

queues = find(~isdelay)';
delays = find(isdelay)';
lGfull = local_lognc(L(queues, :), local_zsum(L, delays, R), N, method);

b = zeros(size(n));
for t = 1:numel(n)
    % numerator: the full constant minus the shells the level set excludes
    idx = find(level <= n(t) - 1);
    corr = -Inf(1, numel(idx));
    for k = 1:numel(idx)
        corr(k) = lGlow(idx(k)) + local_complnc(L, compl, isdelay, N - Mvec(idx(k), :), R, method);
    end
    num = lGfull + log1p(-exp(local_lse(corr) - lGfull));
    % denominator: the flow out of the shell |m| = n-1
    idx = find(level == n(t) - 1);
    den = -Inf;
    for k = 1:numel(idx)
        terms = -Inf(1, R);
        for r = 1:R
            if A(r) <= 0
                continue
            end
            left = N - Mvec(idx(k), :);
            left(r) = left(r) - 1;
            if any(left < 0)
                continue
            end
            terms(r) = log(A(r)) + local_complnc(L, compl, isdelay, left, R, method);
        end
        den = local_lse([den, lGlow(idx(k)) + local_lse(terms)]);
    end
    b(t) = exp(num - den);
end
end

function Z = local_zsum(L, delays, R)
% The aggregate infinite-server demand CLW takes as its think time.
if isempty(delays)
    Z = zeros(1, R);
else
    Z = sum(L(delays, :), 1);
end
end

function lG = local_complnc(L, compl, isdelay, k, R, method)
% log H(k), the complement's constant at a population near the full one.
if any(k < 0)
    lG = -Inf;
    return
end
queues = compl(~isdelay(compl));
delays = compl(isdelay(compl));
lG = local_lognc(L(queues, :), local_zsum(L, delays, R), k, method);
end

function lG = local_lognc(Lq, Z, k, method)
% log G(k) of a set of queue stations with an aggregate think time. Any method
% returning lG at a population vector serves here; CLW is the default because
% its cost does not grow with the population.
if all(k == 0)
    lG = 0;
    return
end
if isempty(Lq)
    % only infinite servers left: G(k) = prod_r Z_r^k_r / k_r!
    lG = 0;
    for r = 1:numel(k)
        if k(r) == 0
            continue
        end
        if Z(r) <= 0
            lG = -Inf;
            return
        end
        lG = lG + k(r)*log(Z(r)) - factln(k(r));
    end
    return
end
if ~strcmpi(method, 'clw')
    line_error(mfilename, 'Only the clw method is wired here; the point evaluation is a plug-in, so another token needs its own call rather than a silent substitution.');
end
[~, lG] = pfqn_clw(Lq, k, Z);
end

function lg = local_open_coeff(rho, isdelayI, kmax)
% log G_I(0..kmax) of an OPEN subnetwork, convolving the per-node series: the
% coefficients of 1/(1-rho z) are rho^k and those of exp(rho z) are rho^k/k!.
lg = [0, -Inf(1, kmax)];
for i = 1:numel(rho)
    li = zeros(1, kmax+1);
    acc = 0;
    for k = 1:kmax
        acc = acc + log(rho(i));
        if isdelayI(i)
            acc = acc - log(k);
        end
        li(k+1) = acc;
    end
    lgnew = -Inf(1, kmax+1);
    for m = 0:kmax
        lgnew(m+1) = local_lse(lg(m+1:-1:1) + li(1:m+1));
    end
    lg = lgnew;
end
end

function Mvec = local_lattice(bound)
R = numel(bound);
stride = ones(1, R);
for r = 2:R
    stride(r) = stride(r-1) * (bound(r-1) + 1);
end
Lsz = prod(bound + 1);
Mvec = zeros(Lsz, R);
for idx = 1:Lsz
    t = idx - 1;
    for r = 1:R
        Mvec(idx, r) = mod(floor(t / stride(r)), bound(r) + 1);
    end
end
end

function lg = local_lgvec(L, isdelayI, Mvec)
% Low-order coefficients of the subnetwork, the only lattice this routine walks.
Lsz = size(Mvec, 1);
lg = -Inf(Lsz, 1);
for idx = 1:Lsz
    m = Mvec(idx, :);
    lg(idx) = local_station_sum(L, isdelayI, m);
end
end

function v = local_station_sum(L, isdelayI, m)
% G_I(m) by direct enumeration of the splits of m across the nodes, which is
% cheap because m is bounded by n-1 and n is small wherever this routine wins.
nodes = size(L, 1);
R = numel(m);
if nodes == 0
    v = 0;
    if any(m > 0)
        v = -Inf;
    end
    return
end
splits = local_splits(m, R);
v = -Inf;
for s = 1:size(splits, 1)
    head = splits(s, :);
    lterm = local_node_term(L(1, :), isdelayI(1), head);
    if lterm == -Inf
        continue
    end
    rest = local_station_sum(L(2:end, :), isdelayI(2:end), m - head);
    v = local_lse([v, lterm + rest]);
end
end

function out = local_splits(m, R)
% Every 0 <= head <= m, one row per vector.
counts = m + 1;
total = prod(counts);
out = zeros(total, R);
for idx = 1:total
    t = idx - 1;
    for r = 1:R
        out(idx, r) = mod(t, counts(r));
        t = floor(t / counts(r));
    end
end
end

function v = local_node_term(Li, isdelayI, k)
% log X_i(k) for one node: multinomial(|k|;k) prod_r L^k_r at a single server,
% prod_r L^k_r/k_r! at an infinite one.
tot = sum(k);
if tot == 0
    v = 0;
    return
end
v = 0;
if ~isdelayI
    v = factln(tot);
end
for r = 1:numel(k)
    if k(r) == 0
        continue
    end
    if Li(r) <= 0
        v = -Inf;
        return
    end
    v = v - factln(k(r)) + k(r)*log(Li(r));
end
end

function s = local_lse(v)
v = v(:)';
m = max(v);
if isempty(m) || ~isfinite(m)
    if isempty(m)
        s = -Inf;
    else
        s = m;
    end
    return
end
s = m + log(sum(exp(v - m)));
end
