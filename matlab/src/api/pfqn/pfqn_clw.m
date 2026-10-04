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
 % g(K) is the coefficient of prod_j z_j^{K_j}, recovered by nested
 % one-dimensional lattice-Poisson inversions (eq. 2.3) with restrictive
 % static scaling (eqs. 5.41-5.46) and log-domain recovery (eq. 7.1).
 %
 % Both of the paper's accelerations are applied:
 %
 %  - DIMENSION REDUCTION BY DECOMPOSITION (Sec. 3, Sec. 5.4). The factors of
 %    (4.5) induce an interdependence graph on the chains (an edge whenever two
 %    chains share a queue). Removing a subset D disconnects the rest into
 %    components S_i(D), and the inversion dimension is |D| + max_i |S_i(D)|
 %    (eq. 3.2), minimized over D (eq. 3.3). With the D variables held on their
 %    contours the remaining factors split into groups with no variable in
 %    common, so each component is inverted separately and the results
 %    multiplied. The dimension reduction sets the inversion ORDER; the scaling
 %    is otherwise unchanged (Sec. 5.4).
 %  - EULER SUMMATION (Sec. 2.4, eqs. 2.20-2.22). The inner sum of (2.3) is
 %    nearly alternating, so for K_j > n+m it is replaced by the Euler sum of
 %    its first n+m+1 terms, applied twice (once for k >= 0, once for k < 0).
 %    Cost per chain drops from 2 l_j K_j to 2 l_j (n+m+1) evaluations, i.e.
 %    prod_j K_j becomes prod_j min(n+m+1, K_j) (eq. 2.26). The order m is
 %    doubled until the paper's own estimate |E(m,n) - E(m,n+1)| falls below
 %    options.euler_tol, and the exact sum is taken once n+m reaches K_j, so
 %    the acceleration never costs accuracy: a fixed 32-term Euler sum is
 %    already 4e-4 nats off at K_j = 200 on Example 8.2.
 %
 % @fn pfqn_clw(L, N, Z, m, options)
 % @param L  (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
 % @param N  (1 x p) closed-chain population vector K.
 % @param Z  (1 x p) aggregate infinite-server relative intensities rho_{j0}
 %                   (think-time term). Default: zeros(1,p).
 % @param m  (q' x 1) queue multiplicities m_i. Default: ones(q',1).
 % @param options struct with optional fields:
 %          .l       (1 x p) inner lattice parameters l_j (roundoff control),
 %                   indexed by chain. Default by inversion DEPTH: 1 at depth 1,
 %                   2 at depths 2-3, 3 deeper (Section 2.2, page 942).
 %          .gamma   (1 x p) aliasing parameters gamma_j (aliasing ~ 10^-gamma_j),
 %                   indexed by chain. Default by depth: 11, 13, 13, 15, ...
 %          .euler   apply Euler summation where K_j > euler_n + euler_m.
 %                   Default true.
 %          .euler_n number of terms summed exactly before averaging (n in
 %                   eq. 2.22). Default 11.
 %          .euler_m starting order of the Euler averaging (m in eq. 2.22).
 %                   Default 20, so 32 terms per half-sum; doubled on demand.
 %          .euler_tol relative tolerance on |E(m,n) - E(m,n+1)| below which the
 %                   Euler sum is accepted. Default 1e-10.
 %          .euler_maxm largest Euler order reached by doubling. Default 160.
 %                   Beyond it the last estimate is returned rather than paying
 %                   for the exact sum.
 %          .beta    (1 x p) multipliers on the scale parameters alpha_j, the
 %                   manual tuning of page 956 (the paper uses 0.8 <= beta <= 1.2
 %                   on its largest examples). Default ones.
 %          .dimred  apply dimension reduction by decomposition. Default true.
 %          .dimred_maxd  largest |D| examined when minimizing (3.3). Default 4.
 % @return G  Normalization constant g(K). Inf if it overflows double range.
 % @return lG Natural logarithm of g(K) (always finite).
 %
 % Validation: reproduces Table I (Example 8.1, p=1, every K up to 2e7) and
 % Table II rows 1-7 (Example 8.2, p=4) to the 7 printed digits, and Table III
 % (Example 8.3, p=11, which the reduction takes to dimension 2 and which is
 % otherwise unreachable) to the same, its last four rows under the paper's own
 % scale tuning, options.beta(1) = 0.8 (0.95 on the last row; page 956). It
 % matches exact convolution (pfqn_ca) to ~1e-9 wherever the magnitudes stay
 % moderate. NOT reproduced: Table II rows 8-9 and Table IV, where the aliasing
 % residual of the outer inversion exceeds the coefficient being extracted --
 % raising gamma_1 by 7 moves the answer by exactly 7 decades, which is the
 % signature of reading the residual rather than the coefficient. The paper
 % computed those two tables with the model-specific ANALYTIC inner inversion of
 % eq. (5.35) ("no l_j is involved for 2 <= j <= 11"), not with the general
 % algorithm, so this is a limit of the scaling of Section 5 rather than of
 % either acceleration: switching Euler summation off changes those rows in no
 % digit at all.
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

useEuler = true;
eulerN = 11;
eulerM = 20;
eulerTol = 1e-10;
eulerMaxM = 160;
useDimred = true;
dimredMaxD = 4;
if isfield(options, 'euler') && ~isempty(options.euler)
    useEuler = logical(options.euler);
end
if isfield(options, 'euler_tol') && ~isempty(options.euler_tol)
    eulerTol = options.euler_tol;
end
if isfield(options, 'euler_maxm') && ~isempty(options.euler_maxm)
    eulerMaxM = round(options.euler_maxm);
end
if isfield(options, 'euler_n') && ~isempty(options.euler_n)
    eulerN = round(options.euler_n);
end
if isfield(options, 'euler_m') && ~isempty(options.euler_m)
    eulerM = round(options.euler_m);
end
if isfield(options, 'beta') && ~isempty(options.beta)
    beta = options.beta(:).';
else
    beta = ones(1, p);
end
if isfield(options, 'dimred') && ~isempty(options.dimred)
    useDimred = logical(options.dimred);
end
if isfield(options, 'dimred_maxd') && ~isempty(options.dimred_maxd)
    dimredMaxD = round(options.dimred_maxd);
end

% trivial populations
if any(N < 0)
    G = 0; lG = -Inf; return
end
if all(N == 0)
    G = 1; lG = 0; return
end

% dimension reduction (Section 3): D is inverted first, then each connected
% component of the interdependence graph minus D, independently
[ordD, comps] = clw_dimred(L, useDimred, dimredMaxD);
order = [ordD, comps{:}];
d = numel(ordD);

% depth of each chain in the nested inversion: D occupies depths 1..d and every
% component restarts at depth d+1, since components are inverted in parallel
depth = zeros(1, p);
depth(ordD) = 1:d;
for c = 1:numel(comps)
    depth(comps{c}) = d + (1:numel(comps{c}));
end

% lattice/aliasing parameters (Section 2.2, page 942), defaulted by depth
if isfield(options, 'l') && ~isempty(options.l)
    l = options.l(:).';
else
    l = 3 * ones(1, p);
    l(depth == 1) = 1;
    l(depth == 2 | depth == 3) = 2;
end
if isfield(options, 'gamma') && ~isempty(options.gamma)
    gam = options.gamma(:).';
else
    gam = 15 * ones(1, p);
    gam(depth == 1) = 11;
    gam(depth == 2 | depth == 3) = 13;
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
% The loop runs in inversion order, which is what dimension reduction changes
% (Section 5.4); chains of distinct components never deflate one another,
% because they share no queue.
alpha = ones(1, p);
used = zeros(qd, 1);          % sum_{k before j} alpha_k rho_{ki} r_k, per queue
etaMat = double(L ~= 0);      % eta_{ki} = 1 iff rho_{ki} ~= 0 (eq. 5.46)
for t = 1:p
    j = order(t);
    Kj = N(j); lj = l(j);
    if Kj == 0
        continue                          % empty chain: no lattice, and 2*lj*Kj = 0 below
    end
    denom = 1 - used;                     % 1 - sum_{k before j} rhobar_{ki}|z_k|
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
        inner = order(t+1:p);             % chains inverted inside chain j
        for n = 1:numel(es)
            qi = qs(n);
            % N_{ij} = mbar_n - 1 + sum_{k inside j} K_k eta_{k,qi} (eq. 5.43)
            Nn = cummb(n) - 1 + sum(N(inner) .* etaMat(qi, inner));
            if Nn <= 0
                an = 1;
            else
                % in the log domain: the product runs over N_{ij} factors below
                % one, and underflows to zero at a few hundred of them, which
                % would silently set alpha_j = 0 and lG = NaN (Example 8.4 has
                % N_{ij} ~ 3000)
                ll = (1:Nn)';
                an = exp(sum(log((Kj + ll) ./ (Kj + 2 * lj * Kj + ll))) / (2 * lj * Kj));
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
    alpha(j) = beta(j) * aj;              % page 956: manual tuning of alpha_j
    aj = alpha(j);
    used = used + aj * L(:, j) * r(j);    % deflate for subsequent chains
end

% split the factors of (4.5) over D and the components: a queue whose chains
% all lie in D is constant during the component inversions, and every other
% queue has all of its non-D chains inside a single component
inD = false(1, p);
inD(ordD) = true;
compOf = zeros(1, p);
for c = 1:numel(comps)
    compOf(comps{c}) = c;
end
qBucket = zeros(qd, 1);       % 0 = D-only, c = component c
for i = 1:qd
    v = find(L(i, :) ~= 0);
    v = v(~inD(v));
    if ~isempty(v)
        qBucket(i) = compOf(v(1));
    end
end

% context for the recursion (local functions, not nested, to avoid MATLAB
% nested-function workspace sharing across recursive calls)
ctx.N = N;
ctx.l = l;
ctx.r = r;
ctx.p = p;
ctx.d = d;
ctx.ordD = ordD;
ctx.comps = comps;
ctx.arho0 = alpha .* Z;     % 1 x p : alpha_j rho_{j0}
ctx.rhoS = L .* alpha;      % q' x p : alpha_j rho_{ji}
ctx.mrow = m;               % q' x 1
ctx.qD = find(qBucket == 0).';
ctx.qC = cell(1, numel(comps));
for c = 1:numel(comps)
    ctx.qC{c} = find(qBucket == c).';
end

% per-group normalization (Section 2.2, page 944: store only ratios of large
% quantities). The scaling above normalizes the WHOLE generating function, not
% each group of factors separately, and a decomposition multiplies the groups'
% inversions together: Example 8.4 splits into ten groups each carrying
% (1-x)^{-100}, individually of modulus e^165 and jointly e^1650, i.e. Inf.
% Every factor of (4.5) has nonnegative coefficients, so the modulus of a group
% over the contour is maximized at w = r; subtracting that exponent holds each
% group at modulus <= 1 and is a constant that cancels exactly in the recovery.
% On the undecomposed path the single group IS the whole function, so the
% offsets stay zero there and that path is unchanged.
ctx.offD = 0;
ctx.offC = zeros(1, numel(comps));
if d > 0 || numel(comps) > 1
    ctx.offD = clw_group_bound(ctx.qD, ordD, r, ctx);
    for c = 1:numel(comps)
        ctx.offC(c) = clw_group_bound(ctx.qC{c}, comps{c}, r, ctx);
    end
end
ctx.euler = useEuler;
ctx.eulerN = eulerN;
ctx.eulerM = eulerM;
ctx.eulerTol = eulerTol;
ctx.eulerMaxM = eulerMaxM;
ctx.chunk = 2e6;            % vectorization chunk for the innermost sum

% run the nested inversion on the scaled generating function -> gbar(K)
gbar = clw_invert_d(1, zeros(1, p), ctx);

% recovery g(K) = prod alpha0_j^{-1} prod alpha_j^{-K_j} gbar(K) (eq. 7.1),
% with the per-group normalization put back
lG = log(gbar) + ctx.offD + sum(ctx.offC) + sum(ctx.arho0) - sum(N .* log(alpha));
if lG > 709
    G = Inf;
else
    G = exp(lG);
end
end

% ---- outer inversion over the committed variables D (Section 3) ----
function val = clw_invert_d(t, W, ctx)
if t > ctx.d
    % D fixed: the remaining factors have no variable in common, so the
    % coefficient of the inner monomial is the product of the components'
    val = clw_gbar_sub(W, ctx.qD, ctx.ordD, ctx.offD, ctx);
    for c = 1:numel(ctx.comps)
        val = val * clw_invert_c(c, 1, W, ctx);
    end
    return
end
j = ctx.ordD(t);
val = clw_lattice(j, W, ctx, @(Wn) clw_invert_d(t + 1, Wn, ctx), false);
if t == 1
    val = real(val);
end
end

% ---- inversion of one component of the interdependence graph minus D ----
function val = clw_invert_c(c, s, W, ctx)
vars = ctx.comps{c};
if s > numel(vars)
    val = clw_gbar_sub(W, ctx.qC{c}, vars, ctx.offC(c), ctx);
    return
end
j = vars(s);
isLeaf = (s == numel(vars));
if isLeaf
    val = clw_lattice(j, W, ctx, @(Wn) clw_gbar_sub(Wn, ctx.qC{c}, vars, ctx.offC(c), ctx), true);
else
    val = clw_lattice(j, W, ctx, @(Wn) clw_invert_c(c, s + 1, Wn, ctx), false);
end
if ctx.d == 0 && s == 1
    % with D empty every component is an independent subnetwork, so its
    % coefficient is real
    val = real(val);
end
end

% ---- one-dimensional lattice-Poisson inversion (eq. 2.3), scaled ----
% Extracts the coefficient of w_j^{K_j}. fn is evaluated at the contour points;
% when vec is true it is called once with all points as rows of W.
function val = clw_lattice(j, W, ctx, fn, vec)
Kj = ctx.N(j); lj = ctx.l(j); rj = ctx.r(j);
if Kj == 0
    % [w_j^0] Gbar = Gbar(w_j=0): the K=0 lattice is the single point 0, and
    % exp(-arho0_j) there cancels the +arho0_j added back into lG
    W(:, j) = 0;
    if vec
        val = sum(fn(W));
    else
        val = fn(W);
    end
    return
end
acc = 0;
for k1 = 0:lj-1
    ph = exp(-1i * pi * k1 / lj);
    acc = acc + ph * clw_inner(j, Kj, lj, k1, rj, W, ctx, fn, vec);
end
val = acc / (2 * lj * Kj * rj^Kj);
end

% ---- inner sum of (2.3) over the lattice index k, with Euler summation ----
% The sum splits at k = 0 into two nearly alternating series (Section 2.4);
% each is replaced by its Euler sum (eq. 2.22). The order m is doubled until
% the paper's own estimate |E(m,n) - E(m,n+1)| falls under the tolerance, and
% the exact sum is taken once m reaches K_j, so accuracy is never traded away.
function inner = clw_inner(j, Kj, lj, k1, rj, W, ctx, fn, vec)
mCur = ctx.eulerM;
while true
    T = ctx.eulerN + mCur;
    if ~ctx.euler || Kj <= T + 1
        kk = (-Kj):(Kj-1);
        v = clw_eval(j, kk, Kj, lj, k1, rj, W, ctx, fn, vec);
        inner = sum(((-1) .^ kk(:)) .* v);
        return
    end
    kk = (-(T+2)):(T+1);
    v = clw_eval(j, kk, Kj, lj, k1, rj, W, ctx, fn, vec);
    vpos = v(T+3:end);                 % k =  s      , s = 0..T+1
    vneg = v(T+2:-1:1);                % k = -(s+1)  , s = 0..T+1
    w1 = clw_euler_weights(ctx.eulerN, mCur);       % s = 0..T
    w2 = clw_euler_weights(ctx.eulerN + 1, mCur);   % s = 0..T+1
    sg = (-1) .^ (0:T+1)';
    E1 = sum(sg(1:T+1) .* w1(:) .* (vpos(1:T+1) - vneg(1:T+1)));
    E2 = sum(sg .* w2(:) .* (vpos - vneg));
    if abs(E1 - E2) <= ctx.eulerTol * abs(E2) || mCur >= ctx.eulerMaxM
        inner = E2;
        return
    end
    mCur = 2 * mCur;
end
end

% ---- evaluate the inverted function at the lattice points of index kk ----
function v = clw_eval(j, kk, Kj, lj, k1, rj, W, ctx, fn, vec)
theta = pi * (k1 + lj * kk) / (lj * Kj);
wj = rj * exp(1i * theta);
nk = numel(wj);
v = zeros(nk, 1);
if vec
    for a = 1:ctx.chunk:nk
        b = min(a + ctx.chunk - 1, nk);
        Wn = repmat(W, b - a + 1, 1);
        Wn(:, j) = wj(a:b).';
        v(a:b) = fn(Wn);
    end
else
    for t = 1:nk
        Wn = W;
        Wn(:, j) = wj(t);
        v(t) = fn(Wn);
    end
end
end

% ---- Euler weights: E(m,n) = sum_i w_i (-1)^i a_i, from eq. (2.22) ----
% E(m,n) = 2^-m sum_{k=0}^{m} C(m,k) S_{n+k} with S_t = sum_{i<=t} (-1)^i a_i,
% so a_i carries the mass of every partial sum that contains it.
function w = clw_euler_weights(n, mm)
b = zeros(1, mm + 1);
b(1) = 1;
for k = 1:mm
    b(k + 1) = b(k) * (mm - k + 1) / k;    % C(mm,k)
end
b = b / 2^mm;
tail = fliplr(cumsum(fliplr(b)));          % tail(k+1) = 2^-mm sum_{j>=k} C(mm,j)
w = ones(1, n + mm + 1);
for i = n + 1:n + mm
    w(i + 1) = tail(i - n + 1);
end
end

% ---- interdependence graph and the minimizing subset D (eqs. 3.1-3.3) ----
function [ordD, comps] = clw_dimred(L, enabled, maxd)
p = size(L, 2);
ordD = zeros(1, 0);
comps = {1:p};
if ~enabled || p <= 2
    return
end
adj = false(p, p);
for i = 1:size(L, 1)
    v = find(L(i, :) ~= 0);
    adj(v, v) = true;                     % each factor is a clique (eq. 3.1)
end
adj(1:p+1:end) = false;
[bestC, bestMax] = clw_components(adj, false(1, p));
best = bestMax;                           % |D| = 0
bestD = zeros(1, 0);
for dd = 1:min(maxd, p - 1)
    if dd >= best
        break                             % dimension is at least |D|
    end
    if nchoosek(p, dd) > 2e5
        break                             % (3.3) is solved by enumeration only
    end
    C = nchoosek(1:p, dd);
    for t = 1:size(C, 1)
        mask = false(1, p);
        mask(C(t, :)) = true;
        [cc, mx] = clw_components(adj, mask);
        if dd + mx < best
            best = dd + mx;
            bestD = C(t, :);
            bestC = cc;
        end
    end
end
if best < p
    ordD = bestD;
    comps = bestC;
end
end

% ---- connected components of the graph with the nodes in mask removed ----
function [comps, mx] = clw_components(adj, mask)
p = size(adj, 1);
lab = zeros(1, p);
nc = 0;
for s = 1:p
    if mask(s) || lab(s) > 0
        continue
    end
    nc = nc + 1;
    stack = s;
    lab(s) = nc;
    while ~isempty(stack)
        v = stack(end);
        stack(end) = [];
        nb = find(adj(v, :) & ~mask & lab == 0);
        lab(nb) = nc;
        stack = [stack, nb]; %#ok<AGROW>
    end
end
comps = cell(1, nc);
mx = 0;
for c = 1:nc
    comps{c} = find(lab == c);
    mx = max(mx, numel(comps{c}));
end
if nc == 0
    comps = cell(1, 0);
end
end

% ---- scaled generating function Gbar restricted to a group of factors ----
% Gbar_S(w) = exp( sum_{j in chains} alpha_j rho_{j0} (w_j - 1) )
%             / prod_{i in queues} (1 - sum_j alpha_j rho_{ji} w_j)^{m_i}
function g = clw_gbar_sub(W, queues, chains, off, ctx)
if isempty(chains)
    expo = zeros(size(W, 1), 1);
else
    expo = (W(:, chains) - 1) * ctx.arho0(chains).';
end
if isempty(queues)
    logden = zeros(size(W, 1), 1);
else
    A = W * ctx.rhoS(queues, :).';        % n x |queues|
    logden = log(1 - A) * ctx.mrow(queues);   % principal complex log
end
g = exp(expo - logden - off);
end

% ---- maximum of log|Gbar_S| over the contour, attained at w = r ----
function off = clw_group_bound(queues, chains, r, ctx)
off = 0;
if ~isempty(chains)
    off = off + sum(ctx.arho0(chains) .* (r(chains) - 1));
end
if ~isempty(queues)
    pole = 1 - ctx.rhoS(queues, :) * r(:);
    pole(pole < realmin) = realmin;
    off = off - sum(ctx.mrow(queues) .* log(pole));
end
end
