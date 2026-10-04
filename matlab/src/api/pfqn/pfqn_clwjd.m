function [G, lG] = pfqn_clwjd(Z, N, mu, visits, lcut, options)
% [G, LG] = PFQN_CLWJD(Z, N, MU, VISITS, LCUT, OPTIONS)
%
% Normalizing constant of a closed product-form network made of an aggregated
% infinite-server (delay) node and an arbitrary number of LIMITED
% JOINT-DEPENDENT (LJD) stations, obtained by numerically inverting the
% multichain generating function with the lattice-Poisson algorithm of
% Choudhury, Leung and Whitt (J. ACM 42(5):935-970, 1995).
%
% This is the joint-dependent generalization of PFQN_CLWOI, which is the
% special case LCUT = 1. It stands to PFQN_CLWOI as PFQN_CLW_LLD stands to
% PFQN_CLW: a per-station cutoff beyond which the rate stops changing turns an
% infinite series into a rational function of the same denominators.
%
% LIMITED JOINT DEPENDENCE. Station i has a rate mu_i(n) that reads the whole
% per-class occupancy vector n (a joint-dependent scaling, sn.jdscaling), but
% saturates coordinatewise: with a cutoff vector l_i = LCUT(i,:),
%
%   mu_i(n) = c_{i,t},   t = t_i(n) = ( min(n_1,l_{i,1}), ..., min(n_R,l_{i,R}) ),
%
% i.e. past l_{i,r} further class-r jobs no longer change the rate. The
% clipped vector t ranges over the finite box prod_r {0,...,l_{i,r}}. Order
% independence is l_i = 1 (t is then the support indicator); a multiserver
% station with c servers is l_i = c, since min(sum n, c) is a function of the
% clipped vector once every l_{i,r} >= c; a load-independent queue is l_i = 1
% with a constant rate.
%
% GENERATING FUNCTION. The multichain generating function factorizes over the
% stations,
%
%   G(z) = exp( sum_r Z_r z_r ) prod_i F_i(z),   F_i(z) = sum_n Phi_i(n) z^n,
%
% with Phi_i the v-weighted balance function of station i,
%
%   Phi_i(0) = 1,   mu_i(n) Phi_i(n) = sum_{r: n_r>0} v_{i,r} Phi_i(n - e_r).
%
% Splitting the count lattice by clipped region, on which mu_i is constant, and
% writing F_{i,t} for the part of F_i carried by the states with t_i(n) = t,
%
%   ( mu_{i,t} - sum_{r: t_r = l_{i,r}} v_{i,r} z_r ) F_{i,t}(z)
%        = sum_{r: t_r >= 1} v_{i,r} z_r F_{i,t-e_r}(z),     F_{i,0} = 1,
%   F_i(z) = sum_{t in box} F_{i,t}(z).                                     (*)
%
% The two sides differ because removing a class-r job leaves the region only
% where the coordinate is UNsaturated: for t_r < l_{i,r} the region fixes
% n_r = t_r, so n - e_r lands in t - e_r; for t_r = l_{i,r} the region is
% n_r >= l_{i,r}, so n - e_r lands in t (n_r > l) or in t - e_r (n_r = l),
% which is what puts v_{i,r} z_r on the left. Hence the singularities are the
% hyperplanes sum_{r in S} v_{i,r} z_r = mu_{i,t} over the SATURATED sets
% S = {r : t_r = l_{i,r}}: at most 2^R of them per station, however large the
% cutoffs are. Setting l_i = 1 reduces (*) to the support recursion of
% PFQN_CLWOI, and R = 1 reduces it to Bertozzi-McKenna eq. 2.19.
%
% INVERSION. G(N) is the coefficient of prod_r z_r^{N_r}, recovered by R nested
% one-dimensional lattice-Poisson inversions (CLW eq. 2.3) on contours of radius
% r_j = 10^{-gamma_j/(2 l_j N_j)}. The restrictive static scaling of CLW eqs.
% 5.41-5.46 runs on the expanded constraint matrix listing one row per
% (station, saturated set S), carrying the unit-pole intensities
% v_{i,r}/min{mu_{i,t} : saturated set of t is S}: the smallest rate over the
% regions sharing a saturated set is the binding one. Rows dominated by a
% superset of no larger rate are dropped first. Recovery is in the log domain
% (CLW eq. 7.1).
%
% SCOPE. The rate must be constant on each clipped region; this is checked on
% probe states and is an error otherwise. Any rate is admissible with
% LCUT = N (the default), the clipping being vacuous on the reachable lattice,
% at the cost given below. Rates that never saturate have no finite rational
% transform and are exactly the case where a large LCUT is mandatory.
%
% COST. prod_r 2 l_r N_r contour points, each costing O(M R prod_r (LCUT(i,r)+1))
% for the M station transforms, against O(M prod_r (N_r+1)(N_r+2)/2) for the
% convolution of PFQN_NCJD. The inversion is linear rather than quadratic in
% each population, but the per-point region box grows with the cutoff: with
% LCUT = N the box is the whole lattice and the convolution wins outright. The
% inversion pays off exactly when the joint dependence saturates early.
%
% Parameters:
%   Z  - (1 x R) think-time demand vector of the aggregated delay node.
%   N  - (1 x R) closed population vector, finite.
%   mu - cell array {1 x M} of function handles, one per LJD station. Each
%        mu{m}(n) returns the total service rate of station m at the per-class
%        occupancy vector n (1 x R). May be empty for a pure delay network.
%   visits - (M x R) matrix, or {1 x M} cell of (1 x R) vectors, of per-station
%        class visit ratios v_{i,r}. Default: unit visits.
%   lcut - (M x R) matrix of per-station per-class saturation cutoffs
%        l_{i,r} >= 1, or a scalar/row broadcast to every station. Entries are
%        clipped to N_r, which is exact: a rate difference at n_r > N_r can only
%        move coefficients with n_r > N_r, and every term reaching z^N has
%        n_r <= N_r at every station. Default: N (no truncation).
%   options - struct with optional fields:
%          .l     (1 x R) inner lattice parameters l_j (roundoff control).
%          .gamma (1 x R) aliasing parameters gamma_j (aliasing ~ 10^-gamma_j).
%          Defaults follow CLW: l_1=1,g_1=11; l_2=l_3=2,g=13; l_j>=4=3,g=15.
%
% Returns:
%   G  - Normalizing constant G(N). Inf if it overflows the double range.
%   lG - log(G(N)) (always finite when G > 0).
%
% Example (delay + a station whose rate saturates at two jobs per class):
%   murate = @(n) 1 + sum(min(n, 2));
%   G = pfqn_clwjd([1 2], [8 6], {murate}, [], [2 2]);
%
% See also PFQN_CLWOI, PFQN_CLW_LLD, PFQN_NCJD, PFQN_MVAJD, PFQN_NCOI.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6
    options = struct();
end
if nargin < 5
    lcut = [];
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
    line_error(mfilename, 'pfqn_clwjd requires finite (closed) populations.');
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

% saturation cutoffs, broadcast and clipped to the reachable lattice
if isempty(lcut)
    L = repmat(N, max(M, 1), 1);
elseif isscalar(lcut)
    L = lcut * ones(max(M, 1), R);
elseif isvector(lcut)
    L = repmat(lcut(:).', max(M, 1), 1);
else
    L = lcut;
end
if M > 0 && (size(L, 1) ~= M || size(L, 2) ~= R)
    line_error(mfilename, 'lcut must be (M x R), a 1 x R row, or a scalar.');
end
L = min(max(round(L), 1), repmat(max(N, 1), size(L, 1), 1));

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
% function restricted to z_r = 0, which kills every F_{i,t} with t_r >= 1. The
% rate handles keep taking full-length occupancy vectors with a zero there.
keep = find(N > 0);
Rk = numel(keep);
Nk = N(keep);
Zk = Z(keep);
lpar = lpar(keep);
gam = gam(keep);
Lk = L(:, keep);

% Per-station region tables over the clipped box prod_r {0,...,l_{i,r}}, in
% mixed radix so that t - e_r always precedes t.
ntreg = ones(1, max(M, 1));
tstride = cell(1, max(M, 1));
muT = cell(1, max(M, 1));
regDec = cell(1, max(M, 1));    % {i}{tl} = [chain, column of t-e_r] rows
regSat = cell(1, max(M, 1));    % {i}{tl} = chains saturated in t
satMask = cell(1, max(M, 1));   % {i}(tl) = bitmask of the saturated set
for m = 1:M
    rad_m = Lk(m, :) + 1;
    st = [1, cumprod(rad_m(1:end-1))];
    tstride{m} = st;
    nt = prod(rad_m);
    ntreg(m) = nt;
    muT{m} = zeros(1, nt);
    regDec{m} = cell(1, nt);
    regSat{m} = cell(1, nt);
    satMask{m} = zeros(1, nt);
    for tl = 0:nt-1
        t = mod(floor(tl ./ st), rad_m);        % 1 x Rk clipped region
        dec = zeros(0, 2);
        sat = [];
        smask = 0;
        for b = 1:Rk
            if t(b) >= 1
                dec(end+1, :) = [b, tl - st(b) + 1]; %#ok<AGROW>
            end
            if t(b) == Lk(m, b)
                sat(end+1) = b; %#ok<AGROW>
                smask = smask + 2^(b-1);
            end
        end
        regDec{m}{tl+1} = dec;
        regSat{m}{tl+1} = sat;
        satMask{m}(tl+1) = smask;
        if tl == 0
            muT{m}(1) = 1;                       % F_0 = Phi(0) = 1, rate unused
        else
            muT{m}(tl+1) = clwjd_regionrate(mu{m}, t, keep, R, N, Lk(m, :), m);
        end
    end
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

% Binding rate of each saturated set: regions sharing a saturated set S share
% the hyperplane sum_{r in S} v_r z_r = mu, so the smallest rate constrains.
nmask = 2 ^ Rk;
muS = inf(max(M, 1), nmask);
for m = 1:M
    for tl = 1:ntreg(m)
        sm = satMask{m}(tl);
        if sm > 0
            muS(m, sm + 1) = min(muS(m, sm + 1), muT{m}(tl));
        end
    end
end

% Restrictive static scaling (CLW eqs. 5.41-5.46) on the expanded constraint
% matrix, one row per (station, saturated set), with dominated rows dropped:
% set S is implied by a superset S' with mu_{i,S'} <= mu_{i,S}, since then
% v_{i,r}/mu_{i,S'} >= v_{i,r}/mu_{i,S} on all of S.
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
Lt = zeros(max(M, 1) * (nmask - 1), Rk);
row = 0;
for m = 1:M
    for mask = 1:nmask-1
        if ~isfinite(muS(m, mask + 1))
            continue
        end
        dominated = false;
        for mask2 = 1:nmask-1
            if mask2 ~= mask && bitand(mask, mask2) == mask ...
                    && isfinite(muS(m, mask2 + 1)) ...
                    && muS(m, mask2 + 1) <= muS(m, mask + 1) * (1 + 1e-12)
                dominated = true;
                break
            end
        end
        if dominated
            continue
        end
        row = row + 1;
        b = find(bitget(mask, 1:Rk) == 1);
        Lt(row, b) = V(m, b) / muS(m, mask + 1);
    end
end
Lt = Lt(1:max(row, 1), :);
nrow = size(Lt, 1);

alpha = ones(1, Rk);
used = zeros(nrow, 1);
etaMat = double(Lt ~= 0);
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
ctx.muT = muT;
ctx.regDec = regDec;
ctx.regSat = regSat;
ctx.ntreg = ntreg;
ctx.chunk = max(1, floor(2e6 / max(ntreg)));

% run the nested inversion on the scaled generating function -> gbar(N)
gbar = clwjd_invert(1, zeros(1, 0), ctx);

% recovery G(N) = exp(sum alpha_r Z_r) prod alpha_r^{-N_r} gbar(N) (eq. 7.1)
lG = log(gbar) + sum(ctx.arho0) - sum(Nk .* log(alpha));
if lG > 709
    G = Inf;
else
    G = exp(lG);
end
end

% ---- rate of one clipped region, with the constancy check ------------------
function rate = clwjd_regionrate(murate, t, keep, R, N, lrow, m)
% The region {n : t(n) = t} pins every unsaturated coordinate and leaves the
% saturated ones free above the cutoff, so the rate is probed at the region
% representative and at two larger occupancies of the same region.
nrep = zeros(1, R);
nrep(keep) = t;
rate = murate(nrep);
if ~(rate > 0)
    line_error(mfilename, sprintf('Station %d has a non-positive rate on a reachable region.', m));
end
sat = find(t == lrow);
if isempty(sat)
    return
end
for pass = 1:2
    nprobe = nrep;
    if pass == 1
        nprobe(keep(sat)) = N(keep(sat));
    else
        nprobe(keep(sat)) = max(t(sat), floor((t(sat) + N(keep(sat))) / 2));
    end
    if any(nprobe ~= nrep)
        rt = murate(nprobe);
        if abs(rt - rate) > 1e-9 * max(1, abs(rate))
            line_error(mfilename, sprintf(['Station %d has a rate that varies within a clipped ' ...
                'region (mu=%g at the region representative, %g above the cutoff). pfqn_clwjd ' ...
                'requires the rate to saturate at lcut; raise lcut or use pfqn_ncjd.'], ...
                m, rate, rt));
        end
    end
end
end

% ---- one-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled ------
% Extracts the coefficient of w_j^{N_j} from g^{(j)}, recursing on inner chains.
function val = clwjd_invert(j, wfixed, ctx)
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
            fv = clwjd_gbar(W, ctx);
            inner = inner + sum(signs(a:b).' .* fv);
        end
    else
        for t = 1:numel(wj)
            inner = inner + signs(t) * clwjd_invert(j + 1, [wfixed, wj(t)], ctx);
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
function g = clwjd_gbar(W, ctx)
% Gbar(w) = exp( sum_r alpha_r Z_r (w_r - 1) ) prod_i F_i(alpha_r v_{i,r} w_r)
% with F_i given by the clipped-region recursion (*). Summing logs is
% legitimate for the principal complex log because exp(log a + log b) = a b.
expo = (W - 1) * ctx.arho0.';              % n x 1
logF = zeros(size(W, 1), 1);
for i = 1:ctx.M
    X = W .* ctx.vs(i, :);                 % n x Rk scaled contour arguments
    nt = ctx.ntreg(i);
    FT = zeros(size(W, 1), nt);
    FT(:, 1) = 1;                          % empty region: Phi(0) = 1
    for tl = 2:nt
        dec = ctx.regDec{i}{tl};
        sat = ctx.regSat{i}{tl};
        num = zeros(size(W, 1), 1);
        for k = 1:size(dec, 1)
            num = num + X(:, dec(k, 1)) .* FT(:, dec(k, 2));
        end
        den = ctx.muT{i}(tl) * ones(size(W, 1), 1);
        for k = 1:numel(sat)
            den = den - X(:, sat(k));
        end
        FT(:, tl) = num ./ den;
    end
    logF = logF + log(sum(FT, 2));
end
g = exp(expo + logF);
end
