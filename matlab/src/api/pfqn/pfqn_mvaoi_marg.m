function [XN, QN] = pfqn_mvaoi_marg(D, N, isDelay, mu)
% [XN, QN] = PFQN_MVAOI_MARG(D, N, ISDELAY, MU)
%
% Exact marginal load-dependent MVA for a closed product-form network of
% infinite-server (delay) and load-independent (single-server, product-form)
% stations plus ANY number of order-independent (OI) stations. An OI station is
% a class-dependent load-dependent server whose total service rate mu_i(n) is a
% permutation-invariant function of the per-class count vector n. The network is
% product-form and is solved exactly by the load-dependent MVA recursion that
% carries, for EACH OI station i, its joint count-vector marginal distribution
% pM_i(n | k):
%
%   pM_i(n | k) = (1/mu_i(n)) * sum_r X_r(k) * pM_i(n - e_r | k - e_r),  n ~= 0
%   pM_i(0 | k) = 1 - sum_{n ~= 0} pM_i(n | k)
%
% All OI-station marginals share the common per-class throughput X_r(k); the
% recursion is exact per station because in product form
% pM_i(n|k) = Phi_i(n) G_{-i}(k-n)/G(k) with X_r(k) = G(k-e_r)/G(k) and the
% balanced-fairness identity Phi_i(n) = (1/mu_i(n)) sum_r Phi_i(n-e_r).
%
% Because the OI rate is class-dependent, the mean-value response-time formula is
% not exact; instead the per-class throughput X_r(k) is closed at each population
% level by population conservation
%
%   X_r(k) * A_r(k) + sum_i QM_ir(k; X) = k_r,   A_r(k) = sum_{i not OI} R_ir(k)
%
% where QM_ir(k) = sum_n n_r pM_i(n | k) is read off each OI station's exact
% marginal. This yields Q, X matching the exact CTMC / normalizing-constant (NC)
% results for any number of OI stations. This is the marginal-distribution
% counterpart of PFQN_MVAOI (the mean-value CMVA form).
%
% Parameters:
%   D       - (M x R) per-class service demand V_ir/rate_ir at every station.
%             Rows of OI stations are ignored (rate comes from MU).
%   N       - (1 x R) closed population vector, finite.
%   isDelay - (1 x M) logical, true for infinite-server (delay) stations.
%   mu      - (1 x M) cell; mu{i} is a function handle mu_i(n) returning the OI
%             total service rate for the per-class occupancy vector n (1 x R) at
%             OI station i, and [] for non-OI stations.
%
% Returns:
%   XN - (1 x R) per-class throughput X_r = G(N-e_r)/G(N).
%   QN - (M x R) per-class mean queue-length at every station.
%
% Reference:
%   Reiser, Lavenberg (1980). Mean-Value Analysis of Closed Multichain Queuing
%   Networks. JACM 27(2). Load-dependent extension: Bruell, Balbo, Afshari
%   (1984). OI stations: Casale, Comte, Dorsman (2026).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = size(D, 1);
R = numel(N);
N = round(N(:)');
isDelay = logical(isDelay(:)');
isOI = false(1, M);
for i = 1:M
    isOI(i) = ~isempty(mu{i});
end
oi_list = find(isOI);
nOI = numel(oi_list);
if nOI == 0
    line_error(mfilename, 'pfqn_mvaoi_marg requires at least one order-independent station.');
end

% OI rank rate mu_i(n) per OI station via the (1-based) canonical microstate.
muM = cell(1, nOI);
for o = 1:nOI
    svcRateFun = mu{oi_list(o)};
    muM{o} = @(n) oi_rate(svcRateFun, n, R);
end

% Caches keyed by population-vector string. The OI marginals are carried per
% station: pM_keys{o}/pM_vals{o} hold station oi_list(o)'s count-vector marginal.
X_cache = configureDictionary('string','cell');   % X_r(k)
Q_cache = configureDictionary('string','cell');   % Q_ir(k)  (M x R)
pM_keys = cell(1, nOI);
pM_vals = cell(1, nOI);
for o = 1:nOI
    pM_keys{o} = configureDictionary('string','cell');
    pM_vals{o} = configureDictionary('string','cell');
end

zeroKey = vec2key(zeros(1,R));
X_cache{zeroKey} = zeros(1,R);
Q_cache{zeroKey} = zeros(M,R);
for o = 1:nOI
    pM_keys{o}{zeroKey} = zeros(1,R);
    pM_vals{o}{zeroKey} = 1.0;
end

pops = enum_vecs(N);
[~, ord] = sortrows([sum(pops,2), pops]);
pops = pops(ord, :);

for pidx = 1:size(pops,1)
    k = pops(pidx, :);
    if sum(k) == 0
        continue
    end
    kkey = vec2key(k);

    % Fixed (non-OI) station response times via the arrival theorem at k - e_r.
    Rfix = zeros(M, R);
    A = zeros(1, R);
    for r = 1:R
        if k(r) == 0
            continue
        end
        kr = k; kr(r) = kr(r) - 1;
        Qkr = Q_cache{vec2key(kr)};
        for i = 1:M
            if isOI(i)
                continue
            end
            if isDelay(i)
                Rfix(i,r) = D(i,r);
            else
                Rfix(i,r) = D(i,r) * (1 + sum(Qkr(i,:)));
            end
            A(r) = A(r) + Rfix(i,r);
        end
    end

    % Close X(k) by population conservation: X_r*A_r + sum_i QM_ir(k;X) = k_r.
    Xk = zeros(1,R);
    for r = 1:R
        if k(r) > 0
            Xk(r) = k(r) / (A(r) + 1);
        end
    end
    nkAll = cell(1, nOI); pvAll = cell(1, nOI);
    for it = 1:2000
        QMtot = zeros(1, R);
        for o = 1:nOI
            [nkAll{o}, pvAll{o}] = oi_marginal(k, Xk, muM{o}, R, pM_keys{o}, pM_vals{o});
            QMtot = QMtot + sum(nkAll{o} .* pvAll{o}, 1);
        end
        Xnew = zeros(1,R);
        for r = 1:R
            if k(r) > 0 && A(r) > 0
                Xnew(r) = max((k(r) - QMtot(r)) / A(r), 0);
            end
        end
        if max(abs(Xnew - Xk)) < 1e-13
            Xk = Xnew;
            break
        end
        Xk = 0.5*Xk + 0.5*Xnew;
    end

    for o = 1:nOI
        [nkAll{o}, pvAll{o}] = oi_marginal(k, Xk, muM{o}, R, pM_keys{o}, pM_vals{o});
    end

    Qk = zeros(M, R);
    for o = 1:nOI
        QMi = sum(nkAll{o} .* pvAll{o}, 1);
        Qk(oi_list(o), :) = QMi;
    end
    for r = 1:R
        for i = 1:M
            if isOI(i)
                continue
            end
            Qk(i,r) = Xk(r) * Rfix(i,r);
        end
    end
    X_cache{kkey} = Xk;
    Q_cache{kkey} = Qk;
    for o = 1:nOI
        pM_keys{o}{kkey} = nkAll{o};
        pM_vals{o}{kkey} = pvAll{o};
    end
end

Nkey = vec2key(N);
XN = X_cache{Nkey};
QN = Q_cache{Nkey};
end

% =========================================================================
% Helper functions
% =========================================================================

function rate = oi_rate(svcRateFun, n, R)
% OI total service rate at count vector n via the 1-based canonical microstate.
if sum(n) == 0
    rate = 0;
    return
end
micro = repelem(1:R, n);
rate = svcRateFun(micro);
end

function [nk, pv] = oi_marginal(k, Xk, muM, R, pM_keys, pM_vals)
% Vector marginal pM(n | k) at the OI station given throughput Xk. Returns count
% vectors nk (rows) and their probabilities pv (column), including n = 0.
vecs = enum_vecs(k);
nrows = size(vecs, 1);
nk = zeros(nrows, R);
pv = zeros(nrows, 1);
cnt = 0;
psum = 0;
zeroRow = -1;
for idx = 1:nrows
    n = vecs(idx, :);
    if sum(n) == 0
        cnt = cnt + 1;
        nk(cnt, :) = n;
        pv(cnt) = 0;         % filled after normalization
        zeroRow = cnt;
        continue
    end
    rate = muM(n);
    if rate <= 0
        continue
    end
    acc = 0;
    for r = 1:R
        if n(r) >= 1 && k(r) >= 1
            nr = n; nr(r) = nr(r) - 1;
            kr = k; kr(r) = kr(r) - 1;
            acc = acc + Xk(r) * lookup_prob(pM_keys, pM_vals, kr, nr, R);
        end
    end
    cnt = cnt + 1;
    nk(cnt, :) = n;
    pv(cnt) = acc / rate;
    psum = psum + pv(cnt);
end
nk = nk(1:cnt, :);
pv = pv(1:cnt);
if zeroRow > 0
    pv(zeroRow) = 1 - psum;
end
end

function p = lookup_prob(pM_keys, pM_vals, kr, nr, R)
% Probability pM(nr | kr) from the cached marginal at population kr.
krkey = vec2key(kr);
if ~isKey(pM_keys, krkey)
    p = 0;
    return
end
keysMat = pM_keys{krkey};
valsVec = pM_vals{krkey};
match = all(keysMat == repmat(nr, size(keysMat,1), 1), 2);
row = find(match, 1);
if isempty(row)
    p = 0;
else
    p = valsVec(row);
end
end

function vecs = enum_vecs(bound)
% All integer vectors v with 0 <= v <= bound (rows), in ndgrid order.
R = numel(bound);
if R == 0
    vecs = zeros(1,0);
    return
end
ranges = cell(1, R);
for r = 1:R
    ranges{r} = 0:bound(r);
end
grids = cell(1, R);
[grids{:}] = ndgrid(ranges{:});
vecs = zeros(numel(grids{1}), R);
for r = 1:R
    vecs(:, r) = grids{r}(:);
end
end

function key = vec2key(v)
key = sprintf('%d_', round(v));
end
