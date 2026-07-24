%{
%{
 % @file pfqn_clw.m
 % @brief Choudhury-Leung-Whitt normalization constant by numerical inversion
 %        of the generating function (JACM 42(5):935-970, 1995).
%}
%}

%{
%{
 % @brief Computes the normalization constant g(K) of a multichain closed
 %        product-form network with single-server and (optionally)
 %        infinite-server queues by numerically inverting its p-dimensional
 %        generating function (Choudhury, Leung and Whitt, 1995).
 %
 % The generating function of g(K) is (eq. 4.5)
 %
 %     G(z) = exp( sum_j rho_{j0} z_j ) / prod_i ( 1 - sum_j rho_{ji} z_j )^{m_i}
 %
 % where j = 1..p indexes closed chains, i = 1..q' indexes the distinct
 % single-server queues with multiplicity m_i, rho_{ji} is the relative
 % traffic intensity of chain j at queue i and rho_{j0} the aggregate
 % relative traffic intensity of chain j at the infinite-server queues.
 % g(K) is the coefficient of prod_j z_j^{K_j}, recovered by p nested
 % one-dimensional lattice-Poisson inversions (eq. 2.3) with restrictive
 % static scaling (eqs. 5.41-5.46) and log-domain recovery (eq. 7.1).
 %
 % @fn pfqn_clw(L, N, Z, m, options)
 % @param L  (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
 % @param N  (1 x p) closed-chain population vector K.
 % @param Z  (1 x p) aggregate infinite-server relative intensities rho_{j0}
 %                   (think-time term). Default: zeros(1,p).
 % @param m  (q' x 1) queue multiplicities m_i. Default: ones(q',1).
 % @param options struct with optional fields:
 %          .l     (1 x p) inner lattice parameters l_j (roundoff control).
 %          .gamma (1 x p) aliasing parameters gamma_j (aliasing ~ 10^-gamma_j).
 %          Defaults follow the paper: l_1=1,g_1=11; l_2=l_3=2,g=13;
 %          l_j>=4 = 3, g=15.
 % @return G  Normalization constant g(K). Inf if it overflows double range.
 % @return lG Natural logarithm of g(K) (always finite).
 %
 % Validation: reproduces the paper's Table I (Example 8.1, p=1, all K up to
 % 2e7) and the direct-summation rows of Table II (Example 8.2, p=4) to 7
 % significant digits, and matches exact convolution (pfqn_ca) to ~1e-9.
 %
 % Scope: this routine implements the exact nested lattice-Poisson inversion
 % with restrictive static scaling. Its cost is prod_j 2 l_j K_j, so it is
 % practical for moderate populations and few chains. The paper's speed-ups
 % for large populations/chains (Euler summation of the inner sums, Section
 % 2.4, and dimension reduction, Sections 3 and 5.4) are not applied here;
 % Tables III-VI and the largest Table II rows rely on them for tractability.
%}
%}
function [G, lG] = pfqn_clw(L, N, Z, m, options)
[qd, p] = size(L);
if nargin < 3 || isempty(Z)
    Z = zeros(1, p);
end
if nargin < 4 || isempty(m)
    m = ones(qd, 1);
end
if nargin < 5
    options = struct();
end
N = N(:).';
Z = Z(:).';
m = m(:);

% default lattice/aliasing parameters (Section 2.2, page 962)
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

% trivial populations
if any(N < 0)
    G = 0; lG = -Inf; return
end
if all(N == 0)
    G = 1; lG = 0; return
end

% contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (eq. 2.7)
r = zeros(1, p);
for j = 1:p
    if N(j) == 0
        r(j) = 1;
    else
        r(j) = 10^(-gam(j) / (2 * l(j) * N(j)));
    end
end

% restrictive static scaling: alpha_j, alpha0_j (eqs. 5.41-5.46), with the
% outer contour variables evaluated at |z_k| = r_k (most restrictive point).
alpha = ones(1, p);
alpha0 = ones(1, p);
used = zeros(qd, 1);          % sum_{k<j} alpha_k rho_{ki} r_k, per queue
etaMat = double(L ~= 0);      % eta_{ki} = 1 iff rho_{ki} ~= 0 (eq. 5.46)
for j = 1:p
    Kj = N(j); lj = l(j);
    denom = 1 - used;                     % 1 - sum_{k<j} rhobar_{ki}|z_k|
    denom(denom <= 0) = eps;
    e = L(:, j) ./ denom;                 % effective intensity per queue
    posq = find(L(:, j) > 0);
    aj = Inf;
    if ~isempty(posq)
        [es, ord] = sort(e(posq), 'descend');
        qs = posq(ord);                   % original queue indices, sorted
        ms = m(qs);                       % aligned multiplicities (mtilde)
        cumrho = cumsum(es) ./ (1:numel(es))';   % rhobar_n (eq. 5.44)
        cummb = cumsum(ms);               % mbar_n (eqs. 5.40/5.45)
        for n = 1:numel(es)
            qi = qs(n);
            % N_{ij} = mbar_n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43)
            Nn = cummb(n) - 1 + sum(N(j+1:p) .* etaMat(qi, j+1:p));
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
    used = used + aj * L(:, j) * r(j);    % deflate for subsequent chains
end

% context for the recursion (local functions, not nested, to avoid MATLAB
% nested-function workspace sharing across recursive calls)
ctx.N = N;
ctx.l = l;
ctx.r = r;
ctx.p = p;
ctx.arho0 = alpha .* Z;     % 1 x p : alpha_j rho_{j0}
ctx.rhoS = L .* alpha;      % q' x p : alpha_j rho_{ji}
ctx.mrow = m;               % q' x 1
ctx.chunk = 2e6;            % vectorization chunk for the innermost sum

% run the nested inversion on the scaled generating function -> gbar(K)
gbar = clw_invert(1, zeros(1, 0), ctx);

% recovery g(K) = prod alpha0_j^{-1} prod alpha_j^{-K_j} gbar(K) (eq. 7.1)
lG = log(gbar) + sum(ctx.arho0) - sum(N .* log(alpha));
if lG > 709
    G = Inf;
else
    G = exp(lG);
end
end

% ---- one-dimensional lattice-Poisson inversion (eq. 2.3), scaled ----
% Extracts the coefficient of w_j^{K_j} from g^{(j)}, recursing on inner chains.
function val = clw_invert(j, wfixed, ctx)
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
            fv = clw_gbar(W, ctx);
            inner = inner + sum(signs(a:b).' .* fv);
        end
    else
        inner = 0;
        for t = 1:numel(wj)
            inner = inner + signs(t) * clw_invert(j + 1, [wfixed, wj(t)], ctx);
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
function g = clw_gbar(W, ctx)
% Gbar(w) = exp( sum_j alpha_j rho_{j0} (w_j - 1) )
%           / prod_i (1 - sum_j alpha_j rho_{ji} w_j)^{m_i}
expo = (W - 1) * ctx.arho0.';             % n x 1
A = W * ctx.rhoS.';                       % n x q'
logden = log(1 - A) * ctx.mrow;           % n x 1 (principal complex log)
g = exp(expo - logden);
end
