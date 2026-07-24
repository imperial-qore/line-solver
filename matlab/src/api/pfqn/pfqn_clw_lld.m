%{
%{
 % @file pfqn_clw_lld.m
 % @brief Choudhury-Leung-Whitt normalization constant by numerical inversion
 %        of the generating function, extended to limited load-dependent (LLD)
 %        stations via the per-center transforms of Bertozzi and McKenna
 %        (SIAM Review 35(2):239-268, 1993).
%}
%}

%{
%{
 % @brief Computes the normalization constant g(K) of a multichain closed
 %        product-form network with limited load-dependent (LLD) stations and
 %        (optionally) infinite-server delay by numerically inverting its
 %        p-dimensional generating function.
 %
 % The generating function is (Bertozzi-McKenna eqs. 2.17/2.23)
 %
 %     G(z) = exp( sum_j rho_{j0} z_j ) prod_i F_i( sum_j rho_{ji} z_j )
 %
 % where F_i is the transform of the station factor of queue i (eq. 2.16),
 %
 %     F_i(x) = sum_{n>=0} x^n / prod_{k=1}^n S_i(k),
 %
 % and S_i(k) = mu(i,k) is the load-dependent rate scaling with k jobs at
 % queue i. For an LLD queue, S_i(k) = c_i constant for k >= l_i, and F_i
 % is the rational function (eq. 2.19)
 %
 %     F_i(x) = [ c_i + sum_{n=1}^{l_i-1} (c_i - S_i(n))
 %                      / prod_{k=1}^n S_i(k) * x^n ] / (c_i - x),
 %
 % analytic except for a simple pole at x = c_i. Multiserver (c_i = number
 % of servers) and load-independent (F_i = 1/(1-x)) queues are special
 % cases. Since g(K) depends on S_i(k) only for k <= sum(K), any general
 % load-dependent input is truncated to LLD at l_i <= sum(K) without loss
 % of exactness.
 %
 % g(K) is the coefficient of prod_j z_j^{K_j}, recovered by p nested
 % one-dimensional lattice-Poisson inversions (CLW, JACM 42(5):935-970,
 % 1995, eq. 2.3) with a restrictive static scaling adapted from CLW eqs.
 % 5.41-5.46: each queue is normalized by its pole c_i (unit-pole form,
 % simple pole) and log-domain recovery (eq. 7.1) is applied.
 %
 % @fn pfqn_clw_lld(L, N, Z, mu, options)
 % @param L  (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
 % @param N  (1 x p) closed-chain population vector K.
 % @param Z  (1 x p) aggregate infinite-server relative intensities rho_{j0}
 %                   (think-time term). Default: zeros(1,p).
 % @param mu (q' x n) load-dependent rate scalings, mu(i,k) = S_i(k); if
 %                    fewer than sum(N) columns are given the last column is
 %                    extended (LLD assumption). Default: ones (all queues
 %                    load-independent).
 % @param options struct with optional fields:
 %          .l     (1 x p) inner lattice parameters l_j (roundoff control).
 %          .gamma (1 x p) aliasing parameters gamma_j (aliasing ~ 10^-gamma_j).
 %          Defaults follow CLW: l_1=1,g_1=11; l_2=l_3=2,g=13; l_j>=4=3,g=15.
 % @return G  Normalization constant g(K). Inf if it overflows double range.
 % @return lG Natural logarithm of g(K) (always finite).
 %
 % Scope: cost is prod_j 2 l_j K_j contour points, each requiring O(sum_i l_i)
 % work, so the routine is practical for moderate populations and few chains.
 % The numerator polynomials are evaluated in double precision; extreme LLD
 % cutoffs (l_i > ~170 with large c_i) may overflow.
%}
%}
function [G, lG] = pfqn_clw_lld(L, N, Z, mu, options)
[qd, p] = size(L);
if nargin < 3 || isempty(Z)
    Z = zeros(1, p);
end
N = N(:).';
Z = Z(:).';
Ntot = sum(N);
if nargin < 4 || isempty(mu)
    mu = ones(qd, max(Ntot, 1));
end
if nargin < 5
    options = struct();
end

% trivial populations
if any(N < 0)
    G = 0; lG = -Inf; return
end
if all(N == 0)
    G = 1; lG = 0; return
end

% extend/truncate mu to sum(N) columns (LLD extension of last column)
if size(mu, 2) < Ntot
    mu = [mu, repmat(mu(:, end), 1, Ntot - size(mu, 2))];
else
    mu = mu(:, 1:Ntot);
end
if any(mu(:) <= 0)
    line_error(mfilename, 'Load-dependent rates mu(i,k) must be positive.');
end

% default lattice/aliasing parameters (CLW Section 2.2, page 962)
if isfield(options, 'l') && ~isempty(options.l)
    l = options.l(:).';
else
    l = 3 * ones(1, p);
    l(1) = 1;
    if p >= 2, l(2) = 2; end
    if p >= 3, l(3) = 2; end
end
if isfield(options, 'gamma') && ~isempty(options.gamma)
    gam = options.gamma(:).';
else
    gam = 15 * ones(1, p);
    gam(1) = 11;
    if p >= 2, gam(2) = 13; end
    if p >= 3, gam(3) = 13; end
end

% drop zero-population chains: the coefficient of z_j^0 equals the pgf
% restricted to z_j = 0, so chain j is removed exactly
keep = N > 0;
L = L(:, keep);
N = N(keep);
Z = Z(keep);
l = l(keep);
gam = gam(keep);
p = numel(N);

% pole c_i and LLD cutoff l_i of each queue: S_i(k) = c_i for k >= l_i
cpole = mu(:, end);
numc = cell(qd, 1);            % numerator coefficients [a_0 ... a_{l_i-1}]
for i = 1:qd
    last = find(mu(i, :) ~= cpole(i), 1, 'last');
    if isempty(last)
        li = 1;                % load-independent up to a constant rate c_i
    else
        li = last + 1;
    end
    a = zeros(1, li);
    a(1) = cpole(i);
    if li > 1
        cp = cumprod(mu(i, 1:li-1));   % prod_{k=1}^n S_i(k)
        a(2:li) = (cpole(i) - mu(i, 1:li-1)) ./ cp;
    end
    numc{i} = a;
end

% contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (CLW eq. 2.7)
r = 10 .^ (-gam ./ (2 * l .* N));

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
Lt = L ./ cpole;               % q' x p, unit-pole intensities
alpha = ones(1, p);
alpha0 = ones(1, p);
used = zeros(qd, 1);           % sum_{k<j} alpha_k rhotilde_{ki} r_k, per queue
etaMat = double(L ~= 0);       % eta_{ki} = 1 iff rho_{ki} ~= 0 (eq. 5.46)
for j = 1:p
    Kj = N(j); lj = l(j);
    denom = 1 - used;                     % 1 - sum_{k<j} rhobar_{ki}|z_k|
    denom(denom <= 0) = eps;
    e = Lt(:, j) ./ denom;                % effective intensity per queue
    posq = find(Lt(:, j) > 0);
    aj = Inf;
    if ~isempty(posq)
        [es, ord] = sort(e(posq), 'descend');
        qs = posq(ord);                   % original queue indices, sorted
        cumrho = cumsum(es) ./ (1:numel(es))';   % rhobar_n (eq. 5.44)
        for n = 1:numel(es)
            qi = qs(n);
            % N_{ij} = n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43, m_i = 1)
            Nn = n - 1 + sum(N(j+1:p) .* etaMat(qi, j+1:p));
            if Nn <= 0
                an = 1;
            else
                ll = (1:Nn)';
                an = (prod((Kj + ll) ./ (Kj + 2 * lj * Kj + ll)))^(1 / (2 * lj * Kj));
            end
            aj = min(aj, an / cumrho(n));
        end
    end
    if Z(j) > 0
        aj = min(aj, Kj / Z(j));          % IS/Poisson term K_j/rho_{j0}
    end
    if ~isfinite(aj)
        aj = 1;                           % chain with no demand anywhere
    end
    alpha(j) = aj;
    alpha0(j) = exp(-aj * Z(j));
    used = used + aj * Lt(:, j) * r(j);   % deflate for subsequent chains
end

% context for the recursion (local functions, not nested, to avoid MATLAB
% nested-function workspace sharing across recursive calls)
ctx.N = N;
ctx.l = l;
ctx.r = r;
ctx.p = p;
ctx.arho0 = alpha .* Z;     % 1 x p : alpha_j rho_{j0}
ctx.rhoS = L .* alpha;      % q' x p : alpha_j rho_{ji}
ctx.q = qd;
ctx.cpole = cpole;          % q' x 1 : F_i pole locations
ctx.numc = numc;            % q' x 1 cell : F_i numerator coefficients
ctx.chunk = 2e6;            % vectorization chunk for the innermost sum

% run the nested inversion on the scaled generating function -> gbar(K)
gbar = clw_lld_invert(1, zeros(1, 0), ctx);

% recovery g(K) = prod alpha0_j^{-1} prod alpha_j^{-K_j} gbar(K) (eq. 7.1)
lG = log(gbar) + sum(ctx.arho0) - sum(N .* log(alpha));
if lG > 709
    G = Inf;
else
    G = exp(lG);
end
end

% ---- one-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled ----
% Extracts the coefficient of w_j^{K_j} from g^{(j)}, recursing on inner chains.
function val = clw_lld_invert(j, wfixed, ctx)
Kj = ctx.N(j); lj = ctx.l(j); rj = ctx.r(j);
acc = 0;
for k1 = 0:lj-1
    ph = exp(-1i * pi * k1 / lj);
    kk = (-Kj):(Kj-1);
    signs = (-1) .^ kk;
    theta = pi * (k1 + lj * kk) / (lj * Kj);
    wj = rj * exp(1i * theta);            % 1 x 2Kj contour points
    if j == ctx.p
        inner = 0;
        nk = numel(wj);
        for a = 1:ctx.chunk:nk
            b = min(a + ctx.chunk - 1, nk);
            W = [repmat(wfixed, b - a + 1, 1), wj(a:b).'];
            fv = clw_lld_gbar(W, ctx);
            inner = inner + sum(signs(a:b).' .* fv);
        end
    else
        inner = 0;
        for t = 1:numel(wj)
            inner = inner + signs(t) * clw_lld_invert(j + 1, [wfixed, wj(t)], ctx);
        end
    end
    acc = acc + ph * inner;
end
val = acc / (2 * lj * Kj * rj^Kj);
if j == 1
    val = real(val);
end
end

% ---- scaled generating function Gbar evaluated at rows of W (n x p) ----
function g = clw_lld_gbar(W, ctx)
% Gbar(w) = exp( sum_j alpha_j rho_{j0} (w_j - 1) )
%           * prod_i F_i( sum_j alpha_j rho_{ji} w_j )
% with F_i(x) = N_i(x)/(c_i - x) (Bertozzi-McKenna eq. 2.19); F_i(0) = 1.
% Note exp(log a + log b) = a*b for the principal complex log, so branch
% choices in the per-queue logs are immaterial.
expo = (W - 1) * ctx.arho0.';              % n x 1
X = W * ctx.rhoS.';                        % n x q' contour arguments
logF = zeros(size(W, 1), 1);
for i = 1:ctx.q
    xi = X(:, i);
    a = ctx.numc{i};
    num = a(end) * ones(size(xi));         % Horner on N_i(x)
    for k = numel(a)-1:-1:1
        num = num .* xi + a(k);
    end
    logF = logF + log(num) - log(ctx.cpole(i) - xi);
end
g = exp(expo + logF);
end
