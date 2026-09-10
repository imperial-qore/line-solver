function [G, lG] = pfqn_clwoi(Z, N, mu, visits, options)
% [G, LG] = PFQN_CLW_OI(Z, N, MU, VISITS, OPTIONS)
%
% Normalizing constant of a closed product-form network made of an aggregated
% infinite-server (delay) node and an arbitrary number of ORDER-INDEPENDENT
% (OI) / pass-and-swap stations with empty swap graph, obtained by numerically
% inverting the multichain generating function with the lattice-Poisson
% algorithm of Choudhury, Leung and Whitt (J. ACM 42(5):935-970, 1995).
%
% This is the OI counterpart of PFQN_CLW_LLD and the transform counterpart of
% the convolution routine PFQN_NCOI. Both return the same G(N); they differ in
% cost, see below.
%
% GENERATING FUNCTION. The multichain generating function factorizes over the
% stations,
%
%   G(z) = exp( sum_r Z_r z_r ) prod_i F_i(z),   F_i(z) = sum_n Phi_i(n) z^n,
%
% with Phi_i the v-weighted balanced-fairness balance function of OI station i,
%
%   Phi_i(0) = 1,   mu_i(n) Phi_i(n) = sum_{r: n_r>0} v_{i,r} Phi_i(n - e_r).
%
% Unlike a load-dependent station, whose factor collapses to a function of the
% single argument sum_r rho_{ri} z_r, an OI station factor genuinely depends on
% the whole vector z, because mu_i(n) depends on the per-class occupancy n
% through its SUPPORT supp(n) = {r : n_r > 0}. It is nevertheless RATIONAL and
% available in closed form. Splitting the count lattice by support, on which
% mu_i(n) = mu_{i,S} is constant, and writing F_{i,S}(z) for the part of F_i
% carried by the states with supp(n) = S, the balance recursion gives
%
%   ( mu_{i,S} - sum_{r in S} v_{i,r} z_r ) F_{i,S}(z)
%        = sum_{r in S} v_{i,r} z_r F_{i,S minus r}(z),      F_{i,{}} = 1,
%   F_i(z) = sum_{S subseteq {1..R}} F_{i,S}(z),                            (*)
%
% since removing one class-r job from a state of support S lands on support S
% when n_r >= 2 and on S minus r when n_r = 1. The singularities are the |S|
% hyperplanes sum_{r in S} v_{i,r} z_r = mu_{i,S}, one per support, in place of
% the single pole x = c_i of the load-dependent case. For R = 1 and constant
% rate c, (*) returns c/(c - v z), and a load-independent single-server queue
% (mu_{i,S} = 1 for all S) returns 1/(1 - sum_r v_{i,r} z_r), i.e. exactly the
% factor PFQN_CLW inverts.
%
% INVERSION. G(N) is the coefficient of prod_r z_r^{N_r}, recovered by R nested
% one-dimensional lattice-Poisson inversions (CLW eq. 2.3) on contours of radius
% r_j = 10^{-gamma_j/(2 l_j N_j)}. The restrictive static scaling of CLW eqs.
% 5.41-5.46 is applied to the EXPANDED constraint matrix that lists one row per
% (station, nonempty support) pair, with unit-pole intensities
% v_{i,r}/mu_{i,S}: each such row is one singular hyperplane of (*), so the
% scaling keeps the whole contour inside the domain of analyticity exactly as
% the single-pole normalization does for PFQN_CLW_LLD. Recovery is in the log
% domain (CLW eq. 7.1).
%
% SCOPE. The rates must be support-only, mu_i(n) = mu_i(supp(n)), which is the
% defining property of an OI station and what makes (*) a finite rational
% function. Every rate handle is therefore verified EXHAUSTIVELY on the count
% lattice 0 < n <= N before the inversion, at prod_r (N_r+1) evaluations per
% station (below the prod_r 2 l_r N_r contour points spent afterwards), and a
% state whose rate differs from that of its support is an error naming that
% state; a rate that varies inside a support is a general balanced-fairness
% station and must go to PFQN_NCOI. The violation is refused rather than
% warned-and-inverted because the inversion would otherwise return a plausible
% but wrong G(N) with no other symptom. A non-empty swap graph breaks the closure of Phi on the
% count vector altogether and requires the microstate routine PFQN_PAS_NC.
% Load-independent single-server queues need no special casing: pass them as OI
% stations with mu_i(n) = 1 and v_{i,r} = D_{i,r}.
%
% COST. prod_r 2 l_r N_r contour points, each costing O(M R 2^R) for the M
% station transforms, against O(M prod_r (N_r+1)(N_r+2)/2) for the convolution
% of PFQN_NCOI. The inversion is therefore LINEAR rather than quadratic in each
% population and wins on large populations with few chains, while the 2^R
% per-point factor makes it lose as the number of chains grows. Unlike
% PFQN_NCOI it returns G at the single population N: throughputs need the R
% additional inversions at N - e_r.
%
% Parameters:
%   Z  - (1 x R) think-time demand vector of the aggregated delay node,
%        Z(r) = 1/sigma_r for a delay with per-class rate sigma_r.
%   N  - (1 x R) closed population vector, finite.
%   mu - cell array {1 x M} of function handles, one per OI station. Each
%        mu{m}(n) returns the total service rate of station m at the per-class
%        occupancy vector n (1 x R) and must depend on n only through its
%        support. May be empty to model a pure delay network.
%   visits - (M x R) matrix, or {1 x M} cell of (1 x R) vectors, of per-station
%        class visit ratios v_{i,r} weighting the balance recursion. Default:
%        unit visits.
%   options - struct with optional fields:
%          .l     (1 x R) inner lattice parameters l_j (roundoff control).
%          .gamma (1 x R) aliasing parameters gamma_j (aliasing ~ 10^-gamma_j).
%          Defaults follow CLW: l_1=1,g_1=11; l_2=l_3=2,g=13; l_j>=4=3,g=15.
%
% Returns:
%   G  - Normalizing constant G(N). Inf if it overflows the double range.
%   lG - log(G(N)) (always finite when G > 0).
%
% Example (delay + one OI station whose rate is the number of busy servers):
%   murate = @(n) sum(n > 0);
%   G = pfqn_clwoi([1 2], [8 6], {murate});
%
% See also PFQN_CLW_LLD, PFQN_CLW, PFQN_NCOI, PFQN_MVAOI, PFQN_PAS_NC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    options = struct();
end
if nargin < 4
    visits = [];
end
if nargin < 3 || isempty(mu)
    mu = {};
end
if ~iscell(mu)
    mu = {mu};
end

R = numel(N);
if isempty(Z)
    Z = zeros(1, R);
end
if numel(Z) ~= R
    line_error(mfilename, 'Z and N must have the same number of classes.');
end
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_clwoi requires finite (closed) populations.');
end
N = round(N(:)');
Z = Z(:)';
M = numel(mu);

if isempty(visits)
    visits = repmat({ones(1, R)}, 1, M);
elseif ~iscell(visits)
    vmat = visits;
    visits = cell(1, M);
    for m = 1:M
        visits{m} = vmat(m, :);
    end
end

% trivial populations
if any(N < 0)
    G = 0; lG = -Inf; return
end
if all(N == 0)
    G = 1; lG = 0; return
end

% default lattice/aliasing parameters (CLW Section 2.2, page 962)
if isfield(options, 'l') && ~isempty(options.l)
    lpar = options.l(:).';
else
    lpar = 3 * ones(1, R);
    lpar(1) = 1;
    if R >= 2, lpar(2) = 2; end
    if R >= 3, lpar(3) = 2; end
end
if isfield(options, 'gamma') && ~isempty(options.gamma)
    gam = options.gamma(:).';
else
    gam = 15 * ones(1, R);
    gam(1) = 11;
    if R >= 2, gam(2) = 13; end
    if R >= 3, gam(3) = 13; end
end

% Drop zero-population chains: the coefficient of z_r^0 is the generating
% function restricted to z_r = 0, which kills every F_{i,S} with r in S. The
% rate handles keep taking full-length occupancy vectors with a zero there.
keep = find(N > 0);
Rk = numel(keep);
Nk = N(keep);
Zk = Z(keep);
lpar = lpar(keep);
gam = gam(keep);
nmask = 2 ^ Rk;

% Support rate table mu_{i,S}, S ranging over the subsets of the retained
% chains encoded as bitmasks 0..2^Rk-1 (mask 0, the empty support, is unused).
muS = zeros(max(M, 1), nmask);
for m = 1:M
    for mask = 1:nmask-1
        chi = zeros(1, R);
        chi(keep(bitget(mask, 1:Rk) == 1)) = 1;
        muS(m, mask + 1) = clwoi_supportrate(mu{m}, chi, m);
    end
end

% Exhaustive support-only verification of every rate handle on the count
% lattice 0 <= n <= N: the transform (*) is exact only if mu_m is constant on
% each support, and a rate that violates that returns a wrong G with no other
% symptom, so it must be refused rather than inverted.
clwoi_checksupport(mu, muS, keep, Nk, R);

% Per-mask chain lists and the S minus r column indices used by recursion (*).
bits = cell(1, nmask - 1);
subcol = cell(1, nmask - 1);
for mask = 1:nmask-1
    b = find(bitget(mask, 1:Rk) == 1);
    bits{mask} = b;
    subcol{mask} = mask - 2 .^ (b - 1) + 1;
end

% Per-station visit vectors restricted to the retained chains.
V = ones(max(M, 1), Rk);
for m = 1:M
    vm = visits{m};
    if numel(vm) ~= R
        line_error(mfilename, 'Each visit vector must have one entry per class.');
    end
    V(m, :) = vm(keep);
end

% contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (CLW eq. 2.7)
rad = 10 .^ (-gam ./ (2 * lpar .* Nk));

% Restrictive static scaling (CLW eqs. 5.41-5.46) on the expanded constraint
% matrix: one row per (station, nonempty support), holding the unit-pole
% intensities v_{i,r}/mu_{i,S} of the singular hyperplane of that support.
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
% Dominated hyperplanes are dropped first: support S of station i is implied by
% a superset S' with mu_{i,S'} <= mu_{i,S}, since then v_{i,r}/mu_{i,S'} >=
% v_{i,r}/mu_{i,S} on all of S. For load-independent stations this leaves the
% single full-support pole and reproduces the scaling of PFQN_CLW_LLD exactly;
% keeping the slack rows would perturb the group averages of eq. 5.44.
Lt = zeros(M * (nmask - 1), Rk);
row = 0;
for m = 1:M
    for mask = 1:nmask-1
        dominated = false;
        for mask2 = 1:nmask-1
            if mask2 ~= mask && bitand(mask, mask2) == mask ...
                    && muS(m, mask2 + 1) <= muS(m, mask + 1) * (1 + 1e-12)
                dominated = true;
                break
            end
        end
        if dominated
            continue
        end
        row = row + 1;
        b = bits{mask};
        Lt(row, b) = V(m, b) / muS(m, mask + 1);
    end
end
Lt = Lt(1:row, :);
nrow = row;
alpha = ones(1, Rk);
used = zeros(nrow, 1);          % sum_{k<j} alpha_k Lt(:,k) r_k, per constraint
etaMat = double(Lt ~= 0);       % eta_{ki} = 1 iff the constraint involves chain k
for j = 1:Rk
    Kj = Nk(j); lj = lpar(j);
    denom = 1 - used;
    denom(denom <= 0) = eps;
    e = Lt(:, j) ./ denom;
    posq = find(Lt(:, j) > 0);
    aj = Inf;
    if ~isempty(posq)
        [es, ord] = sort(e(posq), 'descend');
        qs = posq(ord);
        cumrho = cumsum(es) ./ (1:numel(es))';       % rhobar_n (eq. 5.44)
        for n = 1:numel(es)
            qi = qs(n);
            % N_{ij} = n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43, m_i = 1)
            Nn = n - 1 + sum(Nk(j+1:Rk) .* etaMat(qi, j+1:Rk));
            if Nn <= 0
                an = 1;
            else
                ll = (1:Nn)';
                % in the log domain: the product runs over N_{ij} factors
                % below one, and underflows to zero at a few hundred of them,
                % which would silently set alpha_j = 0 and lG = NaN
                an = exp(sum(log((Kj + ll) ./ (Kj + 2 * lj * Kj + ll))) / (2 * lj * Kj));
            end
            aj = min(aj, an / cumrho(n));
        end
    end
    if Zk(j) > 0
        aj = min(aj, Kj / Zk(j));                    % IS/Poisson term K_j/rho_{j0}
    end
    if ~isfinite(aj)
        aj = 1;                                      % chain with no demand anywhere
    end
    alpha(j) = aj;
    used = used + aj * Lt(:, j) * rad(j);
end

% context for the recursion (local functions, not nested, to avoid MATLAB
% nested-function workspace sharing across recursive calls)
ctx.N = Nk;
ctx.l = lpar;
ctx.r = rad;
ctx.p = Rk;
ctx.M = M;
ctx.arho0 = alpha .* Zk;        % 1 x Rk : alpha_r Z_r
ctx.vs = V .* alpha;            % M x Rk : alpha_r v_{i,r}
ctx.muS = muS;                  % M x 2^Rk support rate table
ctx.bits = bits;
ctx.subcol = subcol;
ctx.nmask = nmask;
ctx.chunk = max(1, floor(2e6 / nmask));

% run the nested inversion on the scaled generating function -> gbar(N)
gbar = clwoi_invert(1, zeros(1, 0), ctx);

% recovery G(N) = exp(sum alpha_r Z_r) prod alpha_r^{-N_r} gbar(N) (eq. 7.1)
lG = log(gbar) + sum(ctx.arho0) - sum(Nk .* log(alpha));
if lG > 709
    G = Inf;
else
    G = exp(lG);
end
end

% ---- support rate of an OI station, read at the support indicator ----------
function rate = clwoi_supportrate(murate, chi, m)
% chi is the 0/1 indicator of the support, itself a lattice point of that
% support since every retained chain has N_r >= 1. Constancy over the support
% is verified afterwards by clwoi_checksupport.
rate = murate(chi);
if ~(rate > 0)
    line_error(mfilename, sprintf('Station %d has a non-positive rate on a reachable support.', m));
end
end

% ---- exhaustive support-only check over the count lattice ------------------
function clwoi_checksupport(mu, muS, keep, Nk, R)
% Every state 0 < n <= N is compared against the rate of its own support. The
% scan costs prod_r (N_r+1) rate evaluations per station, which is below the
% prod_r 2 l_r N_r contour points the inversion itself spends (l_r >= 1 gives
% 2 l_r N_r >= N_r + 1), so exhaustiveness does not change the cost class.
M = numel(mu);
if M == 0
    return
end
Rk = numel(keep);
L = prod(Nk + 1);
for m = 1:M
    n = zeros(1, R);
    for idx = 1:L-1                     % idx 0 is the empty support, unused
        rem_ = idx;
        mask = 0;
        for j = 1:Rk
            nj = mod(rem_, Nk(j) + 1);
            rem_ = floor(rem_ / (Nk(j) + 1));
            n(keep(j)) = nj;
            if nj > 0
                mask = mask + 2 ^ (j - 1);
            end
        end
        rate = mu{m}(n);
        ref = muS(m, mask + 1);
        if abs(rate - ref) > 1e-9 * max(1, abs(ref))
            chi = zeros(1, R);
            chi(keep(bitget(mask, 1:Rk) == 1)) = 1;
            line_error(mfilename, sprintf(['Station %d has a rate that varies within a support: ' ...
                'mu=%g at n=%s but mu=%g at the indicator %s of the same support. ' ...
                'pfqn_clwoi requires order-independent (support-only) rates, mu(n)=mu(supp(n)); ' ...
                'a rate that varies inside a support is a general balanced-fairness station ' ...
                'and must be solved with pfqn_ncoi.'], ...
                m, rate, mat2str(n), ref, mat2str(chi)));
        end
    end
end
end

% ---- one-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled ------
% Extracts the coefficient of w_j^{N_j} from g^{(j)}, recursing on inner chains.
function val = clwoi_invert(j, wfixed, ctx)
Kj = ctx.N(j); lj = ctx.l(j); rj = ctx.r(j);
acc = 0;
for k1 = 0:lj-1
    ph = exp(-1i * pi * k1 / lj);
    kk = (-Kj):(Kj-1);
    signs = (-1) .^ kk;
    theta = pi * (k1 + lj * kk) / (lj * Kj);
    wj = rj * exp(1i * theta);            % 1 x 2Kj contour points
    inner = 0;
    if j == ctx.p
        nk = numel(wj);
        for a = 1:ctx.chunk:nk
            b = min(a + ctx.chunk - 1, nk);
            W = [repmat(wfixed, b - a + 1, 1), wj(a:b).'];
            fv = clwoi_gbar(W, ctx);
            inner = inner + sum(signs(a:b).' .* fv);
        end
    else
        for t = 1:numel(wj)
            inner = inner + signs(t) * clwoi_invert(j + 1, [wfixed, wj(t)], ctx);
        end
    end
    acc = acc + ph * inner;
end
val = acc / (2 * lj * Kj * rj^Kj);
if j == 1
    val = real(val);
end
end

% ---- scaled generating function Gbar evaluated at rows of W (n x Rk) ------
function g = clwoi_gbar(W, ctx)
% Gbar(w) = exp( sum_r alpha_r Z_r (w_r - 1) ) prod_i F_i(alpha_r v_{i,r} w_r)
% with F_i given by the support recursion (*). Summing logs is legitimate for
% the principal complex log because exp(log a + log b) = a b.
expo = (W - 1) * ctx.arho0.';              % n x 1
logF = zeros(size(W, 1), 1);
for i = 1:ctx.M
    X = W .* ctx.vs(i, :);                 % n x Rk scaled contour arguments
    FS = zeros(size(W, 1), ctx.nmask);
    FS(:, 1) = 1;                          % empty support: Phi(0) = 1
    for mask = 1:ctx.nmask-1
        b = ctx.bits{mask};
        sc = ctx.subcol{mask};
        num = zeros(size(W, 1), 1);
        den = ctx.muS(i, mask + 1) * ones(size(W, 1), 1);
        for t = 1:numel(b)
            xt = X(:, b(t));
            num = num + xt .* FS(:, sc(t));
            den = den - xt;
        end
        FS(:, mask + 1) = num ./ den;
    end
    logF = logF + log(sum(FS, 2));
end
g = exp(expo + logF);
end
