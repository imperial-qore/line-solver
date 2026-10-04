function [G, lG, Gtab] = pfqn_ncoi(Z, N, mu, visits, options)
% [G, LG, GTAB] = PFQN_NCOI(Z, N, MU, OPTIONS)
%
% Normalizing constant for a closed product-form queueing network that
% comprises a single aggregated infinite-server (delay) node and an
% arbitrary number of order-independent (OI) / pass-and-swap stations with
% empty swap graph.
%
% The OI stations are analyzed by the balanced-fairness recursion of
% Bonald & Proutiere (2003), "Insensitive bandwidth sharing in data
% networks", combined with the multichain convolution over stations. For a
% single OI station with rank rate mu(supp(n)) the balance function is
%   Phi(0) = 1,  Phi(n) = (1/mu(n)) * sum_{r: n_r>0} Phi(n - e_r),
% and G(N) is obtained by convolving the per-station balance functions with
% the multinomial delay factor F_Z(n) = prod_r Z_r^{n_r} / n_r!,
%   g_0(n) = F_Z(n),   g_m(n) = sum_{0<=x<=n} Phi_m(x) g_{m-1}(n-x),
%   G(N)   = g_M(N).
%
% This is a MACROSTATE routine: both the balance functions and the
% convolution are tabulated over the count lattice 0 <= n <= N, never over
% orderings. That is legitimate exactly because an OI rate is permutation-
% invariant, so Phi(n) -- itself the sum of the ordered-prefix weights
% prod_p 1/mu(c_1..c_p) over all orderings c of the multiset n -- closes on
% the count vector. With a non-empty swap graph that closure fails and the
% microstate routine PFQN_PAS_NC must be used instead.
%
% COST. With L = prod_r (N_r+1) the balance functions cost O(M R L) and the
% convolutions O(M sum_{n<=N} prod_r (n_r+1)) = O(M prod_r (N_r+1)(N_r+2)/2),
% i.e. the order of a load-dependent Buzen convolution: polynomial in the
% population for a fixed number of classes.
%
% Parameters:
%   Z  - (1 x R) think-time demand vector of the aggregated delay node.
%        Z(r) = 1/sigma_r for a delay with per-class rate sigma_r.
%   N  - (1 x R) closed population vector, finite.
%   mu - cell array {1 x M} of function handles, one per OI station. Each
%        mu{m}(n) returns the total service rate of station m given the
%        per-class occupancy (count) vector n (1 x R). For an order-
%        independent station this rate depends only on the support of n
%        (which classes are present), i.e. mu{m}(n) = sum of the capacities
%        of the servers compatible with the classes present in n. A station
%        state whose rate is non-positive is unreachable and is assigned a
%        zero balance value. May be empty to model a pure delay network.
%   options - solver options (optional, currently unused; accepted for
%             signature parity with the other pfqn_* routines).
%
% Returns:
%   G    - Normalizing constant G(N).
%   lG   - log(G(N)).
%   Gtab - (prod_r(N_r+1) x 1) the WHOLE lattice of normalizing constants,
%          Gtab(1 + sum(n .* strides)) = G(n) for every 0 <= n <= N, with
%          strides = [1, cumprod(N(1:end-1)+1)]. The convolution produces this
%          table anyway, so a caller that needs G at more than one population
%          (throughputs, queue lengths, a fold of further stations) must take
%          this output and index it, NOT re-call the routine per population --
%          the latter costs a needless factor prod_r(N_r+1).
%
% Example (IS + two OI stations, R classes):
%   oirate  = @(n) sum(mu1(any(compat1(:, find(n>0)) ~= 0, 2)));
%   oirate2 = @(n) sum(mu2(any(compat2(:, find(n>0)) ~= 0, 2)));
%   G = pfqn_ncoi(1./sigma, N, {oirate, oirate2});
%
% See also PFQN_PAS_NC, PFQN_OI_FNC, PFQN_OI_INSVC, PFQN_MVAOI.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    options = struct(); %#ok<NASGU>
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
    line_error(mfilename, 'pfqn_ncoi requires finite (closed) populations.');
end

N = round(N(:)');
Z = Z(:)';
M = numel(mu);

% Per-OI-station class visit ratios. Default unit visits reproduce the plain
% rank-rate balance; general visits enter as a per-class weight on the
% balanced-fairness recurrence, Phi^v(n)=(1/mu(n)) sum_r v_r Phi^v(n-e_r), the
% v-weighted balance whose class-r geometric factor prod_r v_r^{n_r} carries
% the OI-station visit ratio (see _kb/06-solver-catalog.md, NC OI analyzer).
if nargin < 4 || isempty(visits)
    visits = repmat({ones(1, R)}, 1, M);
elseif ~iscell(visits)
    vmat = visits; visits = cell(1, M);
    for mm = 1:M, visits{mm} = vmat(mm, :); end
end

% Count lattice 0 <= n <= N, flattened column-major; lin(v) = 1 + sum(v.*strides)
% is affine, so lin(x+y) = lin(x) + lin(y) - 1, which the convolution exploits.
dims = N + 1;
ngrid = prod(dims);
strides = [1, cumprod(dims(1:end-1))];
% mixed radix rather than ind2sub, which rejects a single-class (scalar) dims
counts = zeros(ngrid, R);
lin0 = (0:ngrid-1)';
for r = 1:R
    counts(:, r) = mod(floor(lin0 / strides(r)), dims(r));
end
[~, ord] = sort(sum(counts, 2)); % population-increasing sweep order

% Delay balance function: the multinomial factor F_Z(n). A class with
% population but no delay demand makes the state infeasible (weight zero).
g = zeros(ngrid, 1);
for k = 1:ngrid
    v = counts(k, :);
    logf = 0;
    feas = true;
    for r = 1:R
        if v(r) > 0
            if Z(r) <= 0
                feas = false;
                break
            end
            logf = logf + v(r) * log(Z(r)) - gammaln(v(r) + 1);
        end
    end
    if feas
        g(k) = exp(logf);
    end
end

% Convolve in one OI station at a time.
for m = 1:M
    Phim = oi_nc_balance(counts, ord, strides, mu{m}, R, ngrid, visits{m});
    gnext = zeros(ngrid, 1);
    for kx = 1:ngrid
        px = Phim(kx);
        if px == 0
            continue
        end
        x = counts(kx, :);
        rem = N - x;
        % Linear indices of the sub-box 0 <= y <= rem, built by mixed radix.
        ylin = 1;
        for d = 1:R
            ylin = bsxfun(@plus, ylin(:), strides(d) * (0:rem(d)));
            ylin = ylin(:);
        end
        base = 1 + sum(x .* strides);
        gnext(base + ylin - 1) = gnext(base + ylin - 1) + px * g(ylin);
    end
    g = gnext;
end

% g now holds G(n) for EVERY n on the lattice, not just n = N.
Gtab = g;
G = g(ngrid);
lG = log(G);
end

function Phi = oi_nc_balance(counts, ord, strides, murate, R, ngrid, vis)
% v-weighted balanced-fairness recursion
% Phi(n) = (1/mu(n)) sum_{r: n_r>0} v_r Phi(n-e_r), over the count lattice,
% swept in population-increasing order. vis is the 1xR class visit vector.
if nargin < 7 || isempty(vis), vis = ones(1, R); end
Phi = zeros(ngrid, 1);
for t = 1:ngrid
    k = ord(t);
    n = counts(k, :);
    if all(n == 0)
        Phi(k) = 1;
        continue
    end
    mun = murate(n);
    if ~(mun > 0)
        continue % unreachable station state: zero balance value
    end
    acc = 0;
    for r = 1:R
        if n(r) > 0
            acc = acc + vis(r) * Phi(k - strides(r));
        end
    end
    Phi(k) = acc / mun;
end
end
