%{
%{
 % @file pfqn_busyp_multiclass.m
 % @brief Mean busy period of order n for a subnetwork of a multichain network.
%}
%}

%{
%{
 % @brief Mean busy period of order n for a subnetwork of a multichain network.
 % @fn pfqn_busyp_multiclass(alpha, mu, P, N, subnet, n, gamma, phi, tol, jobclass)
 % @param alpha Relative arrival rates (JxR), one column per chain.
 % @param mu Service rates (JxR), the chain-r rate at node j.
 % @param P Routing matrix (JxJ), or a 1xR cell of per-chain matrices.
 % @param N Population per chain (1xR), Inf entries for an open chain.
 % @param subnet Indexes of the nodes forming the subnetwork.
 % @param n Busy period order(s), 1 <= n <= sum(N).
 % @param gamma External arrival rates (JxR), empty for a closed network.
 % @param phi Load-dependent scaling (JxK), dimensionless; empty = single server.
 % @param tol Relative tolerance of the open-network tail truncation.
 % @param jobclass Chain whose own jobs are counted, or -1 to count every chain.
 % @return b Mean busy period duration(s), same size as n.
 % @return lG Log normalizing constants of the subnetwork over the lattice.
 % @return lH Log normalizing constants of the complement over the lattice.
%}
%}
function [b, lG, lH] = pfqn_busyp_multiclass(alpha, mu, P, N, subnet, n, gamma, phi, tol, jobclass)
% [B,LG,LH] = PFQN_BUSYP_MULTICLASS(ALPHA, MU, P, N, SUBNET, N, GAMMA, PHI, TOL, JOBCLASS)
%
% Multichain generalization of PFQN_BUSYP. The busy period of order n for a set
% of nodes I is the time from the instant a job entering I finds n-1 jobs in it
% up to the next instant when fewer than n remain, counting jobs of EVERY chain.
%
% H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean Value
% Analysis", J. ACM 35(3), 1988, states Theorems 1 and 3 for a single chain and
% notes in Section 5 that they carry over to the whole product-form class. The
% proof uses only that the stationary law is product form and that the busy
% period is Keilson's mean ergodic sojourn time on a level set, neither of which
% is single-chain, so replacing the scalar population by a per-chain vector m
% gives, for a closed network,
%
%             sum_{m : |m| >= n}   G_I(m) H(N-m)
%   b(n,I) = --------------------------------------------------
%             sum_{m : |m| = n-1}  G_I(m) sum_r A_r(I) H(N-m-e_r)
%
% with G_I and H the normalizing constants of the subnetwork and of its
% complement at a population VECTOR, and A_r(I) the chain-r arrival flow into I.
% The denominator is the exact chain-r flow across the cut: a chain-r departure
% from the complement at population k occurs at rate alpha_ir H(k-e_r)/H(k), and
% the H(k) cancels the state weight. At R=1 the inner sum holds the single term
% m=n-1 and H(N-m-e_1)=H(N-n), so the expression collapses to Theorem 1 exactly.
%
% THE OPEN CASE NEEDS NO LATTICE. In an open product-form network the stations
% are independent and the total occupancy of node i is a function of the
% AGGREGATE load sum_r alpha_ir/mu_ir alone: summing the station function over
% the compositions of t collapses the multinomial to (sum_r rho_ir)^t. The open
% multichain problem is therefore the single-chain problem on aggregated
% demands, and this routine reduces it to PFQN_BUSYP rather than repeating it.
%
% PER CLASS. With JOBCLASS = r the level set becomes {m_r >= n}, the jobs of
% chain r alone. Only chain-r arrivals move that level, so the flow sum loses its
% sum over r and the same two lattices serve every class.
%
% A MIXED MODEL keeps the closed lattice with its OPEN dimensions TRUNCATED. The
% closed chains are conserved between the subnetwork and its complement, the open
% ones are not: the complement's open count is free, so its open dimensions are
% summed out (Hbar below) and no e_r shift applies to an open chain, removing one
% job from an unbounded dimension leaving the same sum. The truncation is grown
% until the answer stops moving, and is the only approximation in that branch.
%
% For a per-class query on a CLOSED chain of a mixed model there is an exact
% shortcut: marginalizing the open chains leaves a closed network with the
% demands deflated by 1/(1-rho_i^open). Deflate the RATES and not the visits,
% since A_r is built from the visit ratios and folding the deflation into alpha
% would scale the flow constant too.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 7
    gamma = [];
end
if nargin < 8
    phi = [];
end
if nargin < 10 || isempty(jobclass)
    jobclass = -1;
end
if nargin < 9 || isempty(tol)
    % NOT options.tol: that is the iterative-solver tolerance, three orders too
    % loose for a tail sum. Same constant as pfqn_busyp and its ports.
    tol = 1e-12;
end

[J, R] = size(alpha);
N = N(:)';
if numel(N) ~= R
    line_error(mfilename, 'The population vector must have one entry per chain.');
end
isClosed = all(isfinite(N));
isOpen = all(~isfinite(N));
isMixed = ~isClosed && ~isOpen;

subnet = unique(subnet(:)');
if isempty(subnet)
    line_error(mfilename, 'The subnetwork must be non-empty.');
end
if isClosed && numel(subnet) >= J
    line_error(mfilename, 'In a closed network the subnetwork must be a proper subset of the nodes.');
end
if any(subnet < 1) || any(subnet > J)
    line_error(mfilename, 'The subnetwork indexes are out of range.');
end
compl = setdiff(1:J, subnet);

if isempty(phi)
    phi = ones(J, max(1, sum(N(isfinite(N)))));
end

% demands L(i,r) = alpha(i,r)/mu(i,r), zero where chain r does not visit node i
L = zeros(J, R);
visited = alpha > 0 & mu > 0;
L(visited) = alpha(visited) ./ mu(visited);

% A_r(I): the chain-r rate at which jobs enter the subnetwork from outside it
A = zeros(1, R);
for r = 1:R
    Pr = local_routing(P, r);
    A(r) = sum(sum(alpha(compl, r) .* Pr(compl, subnet)));
    if ~isempty(gamma)
        A(r) = A(r) + sum(gamma(subnet, r));
    end
end
if sum(A) <= 0
    line_error(mfilename, 'No job ever enters the subnetwork, its busy period is undefined.');
end
if jobclass > R
    line_error(mfilename, 'The job class index is out of range.');
end
if jobclass >= 1 && A(jobclass) <= 0
    line_error(mfilename, 'No job of that class ever enters the subnetwork, its busy period is undefined.');
end

if isOpen && ~isMixed
    % Exact reduction to the single-chain routine on aggregated demands. The
    % synthetic problem carries no routing, the whole inflow riding on gamma,
    % because A_r(I) is already known and is a constant of the open case.
    rho = sum(L, 2)';
    Pzero = zeros(J, J);
    gsyn = zeros(1, J);
    if jobclass < 1
        % the total occupancy depends on the AGGREGATE load alone
        scalar = rho;
        gsyn(subnet(1)) = sum(A);
    else
        % the class-r marginal is geometric in sigma_ir = rho_ir/(1-rho_i+rho_ir),
        % NOT in rho_ir: the other classes inflate the queue the class-r jobs sit
        % in. That collapse assumes a load-INDEPENDENT station, phi not factoring
        % per class.
        if any(any(phi(subnet, :) ~= 1))
            line_error(mfilename, 'A per-class busy period of an open subnetwork requires load-independent stations: under a load-dependent scaling the class marginal is no longer geometric and the station needs the pair (n_ir,|n_i|) tracked before convolving.');
        end
        denom = 1 - rho + L(:, jobclass)';
        scalar = zeros(1, J);
        scalar(denom > 0) = L(denom > 0, jobclass)' ./ denom(denom > 0);
        gsyn(subnet(1)) = A(jobclass);
    end
    [b, lG] = pfqn_busyp(scalar, phi, Pzero, Inf, subnet, n, gsyn, tol);
    lH = [];
    return
end

if any(n < 1) || any(n ~= round(n))
    line_error(mfilename, 'The busy period order must be a positive integer.');
end
if ~isMixed
    if jobclass < 1
        bound = sum(N);
    else
        bound = N(jobclass);
    end
    if any(n > bound)
        line_error(mfilename, 'The busy period order must be an integer in 1..sum(N), or in 1..N(r) for the busy period of class r alone.');
    end
    [b, lG, lH] = local_lattice_busyp(L, phi, subnet, compl, N, n, A, jobclass, ~isfinite(N));
    return
end

% A mixed model grows the truncation of its open dimensions until the answer
% stops moving; that truncation is the only approximation in this branch.
trunc = 8 + 2*max(n);
prev = [];
while true
    bounds = N;
    bounds(~isfinite(N)) = trunc;
    [b, lG, lH] = local_lattice_busyp(L, phi, subnet, compl, bounds, n, A, jobclass, ~isfinite(N));
    if ~isempty(prev) && all(abs(b - prev) <= 1e-10 * abs(b))
        return
    end
    prev = b;
    trunc = 2*trunc;
    if trunc > 4096
        line_error(mfilename, 'The mixed busy period did not converge: the truncated open dimension keeps growing, so some station of the subnetwork is nearly saturated.');
    end
end
end

function [b, lG, lH] = local_lattice_busyp(L, phi, subnet, compl, bounds, n, A, jobclass, isOpenChain)
% The lattice evaluation shared by the closed and the mixed branch. BOUNDS holds
% the population of a closed chain and the truncation of an open one;
% ISOPENCHAIN names the dimensions that are NOT conserved, whose complement
% counts are summed out rather than read at N-m.
R = numel(bounds);
[Mvec, stride] = local_lattice(bounds);
lG = local_lgvec(L(subnet, :), phi(subnet, :), Mvec, stride, bounds);
lH = local_lgvec(L(compl, :), phi(compl, :), Mvec, stride, bounds);

% Hbar sums the complement over its unconserved dimensions, so it is indexed by
% the CLOSED components alone; with no open chain it is lH itself.
if ~any(isOpenChain)
    lHbar = lH;
else
    keep = ~isOpenChain;
    lHbar = -Inf(size(lH));
    base = local_index(Mvec .* keep, stride);
    for idx = 1:numel(lH)
        lHbar(base(idx)) = local_lse([lHbar(base(idx)), lH(idx)]);
    end
end

% the level set is |m| for the aggregate busy period and m_r for the class-r one;
% only chain-r arrivals move m_r, so the flow sum then holds the single term r
if jobclass < 1
    level = sum(Mvec, 2)';
    chains = 1:R;
else
    level = Mvec(:, jobclass)';
    chains = jobclass;
end
keep = ~isOpenChain;

b = zeros(size(n));
for t = 1:numel(n)
    nt = n(t);
    idx = find(level >= nt);
    rest = (bounds - Mvec(idx, :)) .* keep;
    num = local_lse(lG(idx) + lHbar(local_index(rest, stride)));
    idx = find(level == nt - 1);
    den = -Inf;
    for k = 1:numel(idx)
        terms = -Inf(1, numel(chains));
        for c = 1:numel(chains)
            r = chains(c);
            if A(r) <= 0
                continue
            end
            left = (bounds - Mvec(idx(k), :)) .* keep;
            if ~isOpenChain(r)
                % a closed chain conserves jobs, so the departing one is removed
                left(r) = left(r) - 1;
                if left(r) < 0
                    continue
                end
            end
            terms(c) = log(A(r)) + lHbar(local_index(left, stride));
        end
        den = local_lse([den, lG(idx(k)) + local_lse(terms)]);
    end
    b(t) = exp(num - den);
end
end

function Pr = local_routing(P, r)
% Per-chain routing: one matrix shared by every chain, or one per chain.
if iscell(P)
    Pr = P{r};
else
    Pr = P;
end
end

function [Mvec, stride] = local_lattice(N)
% Every population vector 0 <= m <= N, in the order the linear index implies.
R = numel(N);
stride = ones(1, R);
for r = 2:R
    stride(r) = stride(r-1) * (N(r-1) + 1);
end
Lsz = prod(N + 1);
Mvec = zeros(Lsz, R);
for idx = 1:Lsz
    t = idx - 1;
    for r = 1:R
        Mvec(idx, r) = mod(floor(t / stride(r)), N(r) + 1);
    end
end
end

function idx = local_index(m, stride)
% Linear index of a population vector, one row per vector.
idx = 1 + m * stride';
end

function lg = local_lgvec(L, phi, Mvec, stride, N)
% Log normalizing constants over the whole lattice of a set of nodes.
%
% The station function is
%   X_i(m) = multinomial(|m|; m) prod_r L(i,r)^m_r / prod_{k=1}^{|m|} phi_i(k),
% which at R=1 is the prod_{k} alpha_i/mu_i(k) of pfqn_busyp and at phi(k)=k the
% infinite-server form prod_r L^m_r/m_r!. A node whose scaling row is all ones
% takes the Buzen recursion, O(R) per lattice point; any other node needs the
% full sub-lattice convolution.
[Lnodes, R] = size(L);
Lsz = size(Mvec, 1);
lg = -Inf(Lsz, 1);
lg(1) = 0;
for i = 1:Lnodes
    cols = min(size(phi, 2), max(1, sum(N)));
    isLI = all(phi(i, 1:cols) == 1);
    if isLI
        lgnew = lg;
        for idx = 1:Lsz
            acc = lgnew(idx);
            for r = 1:R
                if Mvec(idx, r) > 0 && L(i, r) > 0
                    acc = local_lse([acc, log(L(i, r)) + lgnew(idx - stride(r))]);
                end
            end
            lgnew(idx) = acc;
        end
        lg = lgnew;
    else
        lX = local_station(L(i, :), phi(i, :), Mvec);
        lgnew = -Inf(Lsz, 1);
        for a = 1:Lsz
            if lg(a) == -Inf
                continue
            end
            for c = 1:Lsz
                if lX(c) == -Inf
                    continue
                end
                s = Mvec(a, :) + Mvec(c, :);
                if all(s <= N)
                    j = local_index(s, stride);
                    lgnew(j) = local_lse([lgnew(j), lg(a) + lX(c)]);
                end
            end
        end
        lg = lgnew;
    end
end
end

function lX = local_station(Li, phii, Mvec)
% log X_i(m) over the lattice for one node.
Lsz = size(Mvec, 1);
R = size(Mvec, 2);
lX = zeros(Lsz, 1);
for idx = 1:Lsz
    m = Mvec(idx, :);
    tot = sum(m);
    v = factln(tot);
    for r = 1:R
        if m(r) > 0
            if Li(r) <= 0
                v = -Inf;
                break
            end
            v = v - factln(m(r)) + m(r) * log(Li(r));
        end
    end
    if isfinite(v) && tot > 0
        cols = size(phii, 2);
        for k = 1:tot
            v = v - log(phii(min(k, cols)));
        end
    end
    lX(idx) = v;
end
end

function s = local_lse(v)
% log-sum-exp, stable when every entry is -Inf
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
