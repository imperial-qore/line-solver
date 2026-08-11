function [X, Qoi, Qli, Qdelay, Soi] = pfqn_mvaoi(Z, N, mu, Dli, visits, options)
% [X, QOI, QLI, QDELAY, SOI] = PFQN_MVAOI(Z, N, MU, DLI, OPTIONS)
%
% Mean-value analysis of a closed product-form queueing network composed of an
% aggregated infinite-server (delay) node, any number of load-independent (LI)
% single-server product-form queues, and any number of order-independent (OI) /
% pass-and-swap stations with empty swap graph. This is the mean-value
% counterpart of PFQN_NCOI and the marginal-distribution form PFQN_MVAOI_MARG:
% it returns the same exact per-class throughput and queue-lengths but WITHOUT
% computing any normalizing constant or joint marginal, using only mean
% quantities (throughputs, demands, queue-lengths) evaluated on shifted models.
% It is the composition-dependent generalization of the Conditional MVA (CMVA)
% of Casale, "A Note on Stable Flow-Equivalent Aggregation in Closed Networks"
% (QUESTA 2009), whose "third form" (rate depending on the full per-class
% occupancy vector) is realized here, extended to MULTIPLE OI stations by
% carrying one rate-shift vector s_i per OI station i.
%
% Throughout, r and s index job classes; i indexes OI stations; j indexes LI
% queues. For a single OI station and no LI queue the analysis recurs on the
% shift vector s_i (the OI occupancy already committed at the bottom of station
% i), Nn = N - s_i the jobs still to distribute:
%   Q^{(S)}(Nn) = sum_r U_r^{(S)}(Nn) ( e_r + Q^{(S+e_r@i)}(Nn - e_r) ),
%   U_r^{(S)}(Nn) = D_r^{(S)}(Nn) X_r^{(S)}(Nn)   (bottom-job utilization),
% with the class-r OI demand and throughput satisfying
%   D_r^{(S)}(Nn) = (1/mu_i(s_i+e_r)) rho_{i,r}^{(S)}(Nn-e_r),                     Nn_r = 1,
%   D_r^{(S)}(Nn) = [X_r^{(S)}(Nn-e_r)/X_r^{(S+e_r@i)}(Nn-e_r)] D_r^{(S)}(Nn-e_r), Nn_r >= 2,
%   rho_{i,r}^{(S)}(M) = rho_{i,r}^{(S)}(M-e_s) X_s^{(S)}(M)/X_s^{(S+e_r@i)}(M),  rho(0)=1, s ~= r,
% and X_r^{(S)}(Nn) closed by population conservation. With K OI stations the
% shift becomes a K x R matrix S (row i = s_i); each OI station keeps its own
% D^i, rho^i and Q^i recursions driven by the common throughput X^{(S)}(Nn), and
% the conservation identity aggregates every station's contribution:
%   Nn_r = X_r Z_r + sum_j Q^{(j)}_r + sum_i Q^{(i)}_r,
% where the LI queue Q^{(j)}_r = X_r D_{j,r} (1 + sum_s Q^{(j)}_s(Nn - e_r)) is
% the standard arrival-theorem term. States (S, Nn) are processed by increasing
% sum(Nn) so every reference lands at a strictly smaller free population.
%
% Parameters:
%   Z   - (1 x R) think-time demand vector of the aggregated delay node.
%   N   - (1 x R) closed population vector, finite.
%   mu  - cell array {mu_1,...,mu_K} of function handles; mu_i(n) returns the OI
%         total service rate of station i for per-class occupancy n (1 x R). A
%         bare function handle is accepted as the single-station shorthand.
%   Dli - (J x R) per-class demand matrix of the LI single-server queues
%         (D_{j,r} = V_{j,r}/rate_{j,r}); empty or omitted when J = 0.
%   options - solver options (optional, currently unused).
%
% Returns:
%   X      - (1 x R) per-class throughput X_r = G(N-e_r)/G(N).
%   Qoi    - (K x R) per-class mean queue-length at each OI station (row i).
%   Qli    - (J x R) per-class mean queue-length at each LI queue (row j).
%   Qdelay - (1 x R) per-class mean queue-length at the delay node (X.*Z).
%   Soi    - (K x R) per-class mean number of IN-SERVICE jobs at each OI station,
%            i.e. E[sir_r] with sir_r the count of class-r jobs receiving a
%            strictly positive rank rate (see PFQN_OI_INSVC); the utilization of
%            OI station i is Soi(i,r)/c_i. Unlike X/Qoi/Qli, which are pure
%            mean-value quantities, Soi is a distributional statistic and is
%            therefore obtained from the OI count marginal
%              pM_i(n|k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n-e_r|k-e_r),
%              pM_i(0|k) = 1 - sum_{n ~= 0} pM_i(n|k),
%            which is assembled here from the zero-shift throughputs X^{(0)}(k)
%            already cached by the mean-value recursion above (no normalizing
%            constant is formed). It is only computed when requested.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    Dli = [];
end
if nargin < 6
    options = struct(); %#ok<NASGU>
end
if isa(mu, 'function_handle')
    mu = {mu};
end
if ~iscell(mu) || isempty(mu)
    line_error(mfilename, 'mu must be a (nonempty) cell of OI rate handles.');
end
K = numel(mu);
% Per-OI-station class visit ratios v_{i,r}. They enter the class-r demand base
% case theta_{i,r}(N_r=1)=v_{i,r}/mu_i(...) (eq mvaoi-D); the N_r>=2 ratio case
% cancels visits. Default unit visits. ms-promoted OI stations pass ones here,
% their visits already folded into the rate handle by the analyzer.
if nargin < 5 || isempty(visits)
    visits = repmat({ones(1, numel(N))}, 1, K);
end
for i = 1:K
    if ~isa(mu{i}, 'function_handle')
        line_error(mfilename, 'each mu{i} must be a function handle mu_i(n).');
    end
end

R = numel(N);
N = round(N(:)');
Z = Z(:)';
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_mvaoi requires finite (closed) populations.');
end
if isempty(Dli)
    Dli = zeros(0, R);
end
J = size(Dli, 1);
I = eye(R);

% ---- compact dense state index -------------------------------------------
% A state is (S, Nn): shift matrix S (K x R) plus free population Nn (1 x R),
% subject to sum_i S(i,r) + Nn(r) <= N(r). That constraint separates over
% classes, so the class-r column (S(1,r),...,S(K,r),Nn(r)) is an arbitrary
% (K+1)-vector of sum <= N(r). Ranking those per class (the sum <= N(r) slice
% of MULTICHOOSE, i.e. what SPROD enumerates with its slack row) and mixing
% radix over classes gives an O(1) integer index over exactly the reachable
% set, replacing the earlier SPRINTF-keyed CONTAINERS.MAP lookups.
Vt = cell(1, R); Lr = zeros(1, R); brad = cell(1, R); rnk = cell(1, R);
for r = 1:R
    Mtab = multichoose(K+2, N(r));
    Vt{r} = Mtab(:, 1:K+1);             % drop the slack entry: sum <= N(r)
    Lr(r) = size(Vt{r}, 1);
    brad{r} = (N(r)+1).^(0:K);
    rnk{r} = zeros((N(r)+1)^(K+1), 1);
    rnk{r}(1 + Vt{r}*brad{r}(:)) = 1:Lr(r);
end
cstride = ones(1, R);
for r = 2:R
    cstride(r) = cstride(r-1) * Lr(r-1);
end
nst = prod(Lr);

% Per-class rank transitions: decN drops one free class-r job; incT commits one
% class-r job to station i without freeing one; shfT does both at once, i.e.
% (S,Nn) -> (S + e_r@i, Nn - e_r). A zero entry means out of the simplex.
decN = cell(1, R); incT = cell(1, R); shfT = cell(1, R);
for r = 1:R
    decN{r} = zeros(Lr(r), 1);
    incT{r} = zeros(Lr(r), K);
    shfT{r} = zeros(Lr(r), K);
    for a = 1:Lr(r)
        v = Vt{r}(a, :);
        if v(K+1) > 0
            w = v; w(K+1) = w(K+1) - 1;
            decN{r}(a) = rnk{r}(1 + w*brad{r}(:));
        end
        for i = 1:K
            if sum(v) < N(r)
                w = v; w(i) = w(i) + 1;
                incT{r}(a,i) = rnk{r}(1 + w*brad{r}(:));
            end
            if v(K+1) > 0
                w = v; w(i) = w(i) + 1; w(K+1) = w(K+1) - 1;
                shfT{r}(a,i) = rnk{r}(1 + w*brad{r}(:));
            end
        end
    end
end

% Decode every state once: per-class rank and free population.
Arank = zeros(nst, R);
Nnall = zeros(nst, R);
for r = 1:R
    Arank(:,r) = mod(floor((0:nst-1)'/cstride(r)), Lr(r)) + 1;
    Nnall(:,r) = Vt{r}(Arank(:,r), K+1);
end
[~, ord] = sort(sum(Nnall, 2));         % population-increasing sweep order

% Zero-initialised storage. States with sum(Nn) = 0 keep those zeros for EVERY
% shift S, which is what the Qsub reference (S + e_s@i, Nn - e_s) needs when
% Nn = e_s; no explicit base case is required.
Xt  = zeros(nst, R);        % X^{(S)}(Nn)
Qlt = zeros(nst, J*R);      % Qli^{(S)}(Nn), the J x R matrix flattened
Dt3 = zeros(nst, R, K);     % D_i^{(S)}(Nn)
Qt3 = zeros(nst, R, K);     % Q_i^{(S)}(Nn)
Rh3 = nan(nst, R, K);       % rho_i^{(S)}(M), NaN until computed

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
for p = 1:nst
    k = ord(p);
    Nn = Nnall(k, :);
    if sum(Nn) == 0
        continue
    end
    ak = Arank(k, :);
    S = zeros(K, R);
    for r = 1:R
        S(:,r) = Vt{r}(ak(r), 1:K)';
    end
    kdec = zeros(1, R);                 % index of (S, Nn - e_r)
    for r = 1:R
        if Nn(r) > 0
            kdec(r) = k + (decN{r}(ak(r)) - ak(r)) * cstride(r);
        end
    end

    % Per-OI-station class demands D_i and queue-ahead vectors Qsub_i.
    Dt   = zeros(K, R);         % Dt(i,r)
    Qsub = zeros(R, R, K);      % Qsub(s,r,i) = Q_i^{(S+e_s@i)}(Nn - e_s), comp r
    for i = 1:K
        for r = 1:R
            if Nn(r) == 0
                continue
            end
            if Nn(r) == 1
                mur = mu{i}(S(i,:) + I(r,:));
                if mur > 0
                    Dt(i,r) = (visits{i}(r)/mur) * rho_chain(i, r, kdec(r));
                end
            else
                kshf = k + (shfT{r}(ak(r), i) - ak(r)) * cstride(r);
                if Xt(kshf, r) > 0
                    Dt(i,r) = (Xt(kdec(r), r)/Xt(kshf, r)) * Dt3(kdec(r), r, i);
                end
            end
        end
        for s = 1:R
            if Nn(s) > 0
                kshf = k + (shfT{s}(ak(s), i) - ak(s)) * cstride(s);
                Qsub(s,:,i) = Qt3(kshf, :, i);
            end
        end
    end

    % LI-queue arrival-theorem coefficients beta_j,r and aggregate per class.
    betaLI = zeros(J, R);       % beta_j,r = D_j,r (1 + sum_s Qli_j,s(Nn - e_r))
    for r = 1:R
        if Nn(r) == 0 || J == 0
            continue
        end
        Qsub_li = reshape(Qlt(kdec(r), :), J, R);
        betaLI(:,r) = Dli(:,r) .* (1 + sum(Qsub_li, 2));
    end

    % Population conservation  A X = Nn  over classes with Nn_r > 0.
    idx = find(Nn > 0);
    m = numel(idx);
    A = zeros(m, m);
    for a = 1:m
        r = idx(a);
        for b = 1:m
            s = idx(b);
            if s == r
                val = Z(r) + sum(betaLI(:,r));
                for i = 1:K
                    val = val + Dt(i,r) * (1 + Qsub(r,r,i));
                end
                A(a,b) = val;
            else
                val = 0;
                for i = 1:K
                    val = val + Dt(i,s) * Qsub(s,r,i);
                end
                A(a,b) = val;
            end
        end
    end
    % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
    Xk = zeros(1, R);
    for a = 1:m
        r = idx(a);
        denom = A(a,a);
        for b = 1:m
            if b ~= a
                s = idx(b);
                % Xt(kdec(r),s) = X_s(S, Nn-e_r);  Xt(kdec(s),r) = X_r(S, Nn-e_s)
                if Xt(kdec(s), r) > 0
                    denom = denom + A(a,b) * (Xt(kdec(r), s) / Xt(kdec(s), r));
                end
            end
        end
        if denom > 0
            Xk(r) = Nn(r) / denom;
        end
    end

    % Assemble OI-station queues Q_i and LI-queue queues Qli.
    Qk_li = zeros(J, R);
    for r = 1:R
        if Nn(r) == 0
            continue
        end
        Qk_li(:,r) = Xk(r) * betaLI(:,r);
    end
    for i = 1:K
        U = Dt(i,:) .* Xk;               % 1 x R bottom-job utilization
        Qi = zeros(1, R);
        for r = 1:R
            Qi(r) = U(r) + U * Qsub(:,r,i);
        end
        Qt3(k,:,i) = Qi;
        Dt3(k,:,i) = Dt(i,:);
    end
    Xt(k,:)  = Xk;
    Qlt(k,:) = Qk_li(:)';
end

% Full-population state: every class at rank of [0,...,0,N(r)].
kN = 1;
for r = 1:R
    aN = rnk{r}(1 + N(r)*brad{r}(K+1));
    kN = kN + (aN - 1) * cstride(r);
end
X = Xt(kN, :);
Qoi = zeros(K, R);
for i = 1:K
    Qoi(i,:) = Qt3(kN, :, i);
end
Qli = reshape(Qlt(kN,:), J, R);
Qdelay = X .* Z;

if nargout >= 5
    % X^{(0)}(n) over the population lattice 0 <= n <= N, read off the
    % zero-shift slice of the dense index.
    shp = N + 1;
    lstride = ones(1, R);
    for d = 2:R
        lstride(d) = lstride(d-1) * shp(d-1);
    end
    Xlat = zeros(prod(shp), R);
    n = pprod(N);
    while n >= 0
        kz = 1;
        for r = 1:R
            kz = kz + (rnk{r}(1 + n(r)*brad{r}(K+1)) - 1) * cstride(r);
        end
        Xlat(1 + sum(n .* lstride), :) = Xt(kz, :);
        n = pprod(n, N);
    end
    Soi = oi_insvc_means(N, mu, Xlat, K, R);
end

    function rv = rho_chain(i, r, kstart)
    % rho_{i,r}^{(S)}(M) by walking M -> M - e_s (s ~= r, M_s > 0) down to
    % M = 0, then unwinding. Every visited state is memoised in Rh3, so the
    % chain stops early at the first state already computed.
    chain = zeros(1, sum(N) + 1);
    ratios = zeros(1, sum(N) + 1);
    nc = 0;
    kk = kstart;
    while true
        if ~isnan(Rh3(kk, r, i))
            rv = Rh3(kk, r, i);
            break
        end
        if sum(Nnall(kk,:)) == 0
            rv = 1.0;
            Rh3(kk, r, i) = rv;
            break
        end
        s = 0;
        for ss = 1:R
            if ss ~= r && Nnall(kk, ss) > 0
                s = ss;
                break
            end
        end
        akk = Arank(kk, :);
        k2 = kk + (incT{r}(akk(r), i) - akk(r)) * cstride(r);
        ratio = 0.0;
        if Xt(k2, s) > 0
            ratio = Xt(kk, s) / Xt(k2, s);
        end
        nc = nc + 1;
        chain(nc) = kk;
        ratios(nc) = ratio;
        kk = kk + (decN{s}(akk(s)) - akk(s)) * cstride(s);
    end
    for t = nc:-1:1
        rv = rv * ratios(t);
        Rh3(chain(t), r, i) = rv;
    end
    end
end

% =========================================================================
function Soi = oi_insvc_means(N, mu, Xlat, K, R)
% Mean number of in-service jobs per class at each OI station, E[sir_r], from
% the OI count marginal pM_i(n|k) built on the cached zero-shift throughputs
% X^{(0)}(k). Exact because in product form
% pM_i(n|k) = Phi_i(n) G_{-i}(k-n)/G(k) and X_r(k) = G(k-e_r)/G(k), so the
% recursion below reproduces the balanced-fairness identity for Phi_i.
shp = N + 1;
stride = ones(1, R);
for d = 2:R
    stride(d) = stride(d-1) * shp(d-1);
end
total = prod(shp);
subs = zeros(total, R);
for i = 1:total
    li = i - 1;
    for d = 1:R
        subs(i, d) = mod(li, shp(d));
        li = floor(li / shp(d));
    end
end
[~, ord] = sort(sum(subs, 2));      % process populations by increasing size

Xk = Xlat;                          % X^{(0)}(k) over the lattice
Soi = zeros(K, R);
for m = 1:K
    gm = pfqn_oi_insvc(mu{m}, N);        % E[sir_r | n] over the lattice
    muv = zeros(total, 1);
    for i = 1:total
        if sum(subs(i,:)) > 0
            muv(i) = mu{m}(subs(i,:));
        end
    end
    % pMv(a,b) = pM_m(n_a | k_b), filled for n_a <= k_b by increasing sum(k).
    pMv = zeros(total, total);
    pMv(1,1) = 1;                        % pM(0|0) = 1
    for bb = 1:total
        b = ord(bb);
        k = subs(b, :);
        if sum(k) == 0
            continue
        end
        acc0 = 0;
        for aa = 1:total
            a = ord(aa);
            n = subs(a, :);
            if sum(n) == 0 || any(n > k)
                continue
            end
            if muv(a) <= 0
                continue
            end
            acc = 0;
            for r = 1:R
                if n(r) > 0
                    acc = acc + Xk(b, r) * pMv(a - stride(r), b - stride(r));
                end
            end
            pMv(a, b) = acc / muv(a);
            acc0 = acc0 + pMv(a, b);
        end
        pMv(1, b) = 1 - acc0;            % empty-state probability by complement
    end
    idxN = 1 + sum(N .* stride);
    for r = 1:R
        Soi(m, r) = pMv(:, idxN)' * gm(:, r);
    end
end
end

