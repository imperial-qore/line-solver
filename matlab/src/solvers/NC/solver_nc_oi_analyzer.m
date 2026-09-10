function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_oi_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_OI_NC_ANALYZER(SN, OPTIONS)
%
% Exact normalizing-constant analysis of a closed queueing network that mixes
% order-independent (OI) stations with ordinary BCMP product-form stations.
% Supported stations:
%   - OI stations (SchedStrategy.OI / PAS with an empty swap graph), analyzed
%     by the balanced-fairness rank rate mu(supp n) (Bonald & Proutiere 2003);
%   - any BCMP product-form station: infinite-server (delay, IS), processor
%     sharing (PS), LCFS-PR, and class-independent-rate FCFS, single- or multi-
%     server, analyzed by the load-dependent BCMP weight table
%       W_i(n) = (sum n)!/prod(n_r!) * prod_r D_{i,r}^{n_r} / prod_{k=1}^{sum n} beta_i(k),
%     with D_{i,r} = V(i,r)/rate(i,r) the per-class demand and beta_i(k) the
%     load-dependent capacity (min(k,c) for a c-server queue, k for IS).
%
% The full-network normalizing-constant table G(P) over the lattice
% 0 <= P <= N is assembled by balanced-fairness convolution of all station
% tables, with the OI stations and the aggregated delay evaluated through
% PFQN_NCOI. The exact per-class mean queue length at any station follows from
% the OI functional-server (FNC) identity of PFQN_OI_FNC (Casale, QEST 2006):
%   E[n_{i,r}] = ( sum_{0<=b<=N} Psi_{i,r}(b) G(N-b) ) / G(N) - 1,
% Psi_{i,r} being the FNC balance function built from that station's balance
% table with f(n)=n_r. Per-class throughput X_r = G(N-e_r)/G(N); delay queue
% length and response time follow from Little's law.
%
% General per-class visits at the OI stations are supported: they enter the
% v-weighted balanced-fairness balance Phi^v(n)=(1/mu(n)) sum_r v_r Phi^v(n-e_r).
% BCMP-station visits are arbitrary (folded into D).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;
iter = 1;
method = 'oi';

M = sn.nstations;
K = sn.nclasses;

% ---- the route's own premises (class switching, open chains, a rate lattice,
% a fork, a missing rank rate) are decided by NC_OI_REFUSAL, which SolverNC's
% support gate asks too, so a pair the report offers is a pair this analyzer
% accepts ------------------------------------------------------------------
oiReason = nc_oi_refusal(sn, 'oi');
if ~isempty(oiReason)
    line_error(mfilename, oiReason);
end
N = round(sn.njobs(:)');

% ---- classify stations -----------------------------------------------------
% OI: rank-rate balanced-fairness station. INF: aggregated into the delay Z.
% Q : ordinary BCMP product-form station (PS / LCFS-PR / FCFS), a load-
% dependent weight table with server count c(ist).
isOI  = false(M,1);
isINF = false(M,1);
isQ   = false(M,1);
svc   = cell(M,1);
for ist = 1:M
    ind = sn.stationToNode(ist);
    if sn.sched(ist) == SchedStrategy.INF
        isINF(ist) = true;
    elseif sn.sched(ist) == SchedStrategy.PAS || sn.sched(ist) == SchedStrategy.OI
        sg = [];
        if ind >= 1 && ind <= numel(sn.nodeparam) && isstruct(sn.nodeparam{ind}) ...
                && isfield(sn.nodeparam{ind}, 'swapGraph')
            sg = sn.nodeparam{ind}.swapGraph;
        end
        if isempty(sg) || any(sg(:) ~= 0)
            line_error(mfilename, 'solver_nc_oi supports OI stations only (PAS with a non-empty swap graph is not order-independent).');
        end
        isOI(ist) = true;
        svc{ist} = sn.nodeparam{ind}.svcRateFun;
        if isempty(svc{ist})
            line_error(mfilename, 'OI station %d has no service rate function; set it via setService(@(c) ...).', ist);
        end
    elseif any(sn.sched(ist) == [SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.FCFS, SchedStrategy.SIRO])
        isQ(ist) = true;
        if any(sn.sched(ist) == [SchedStrategy.FCFS, SchedStrategy.SIRO])
            % BCMP type 1 (and the order-insensitive SIRO, which shares the
            % FCFS queue-length distribution for exponential service) require a
            % class-independent service rate for product form.
            rr = sn.rates(ist, :);
            rr = rr(isfinite(rr) & sn.njobs > 0);
            if ~isempty(rr) && (max(rr) - min(rr)) > 1e-9 * max(rr)
                line_error(mfilename, 'Station %d has class-dependent FCFS/SIRO rates and is not product form; solver_nc_oi requires class-independent rates.', ist);
            end
        end
    else
        line_error(mfilename, 'solver_nc_oi supports only INF (delay), OI, PS, LCFS-PR, SIRO and class-independent FCFS stations.');
    end
end

% ---- per-class visits (chain == class); normalize to the reference station -
V = zeros(M, K);
for r = 1:K
    c = find(sn.chains(:, r));           % the chain carrying class r
    vis = sn.visits{c};                  % (nstateful x nclasses)
    for ist = 1:M
        isf = sn.stationToStateful(ist);
        V(ist, r) = vis(isf, r);
    end
    vref = V(sn.refstat(r), r);
    if vref > 0
        V(:, r) = V(:, r) / vref;
    end
end

% ---- per-class demand and aggregated delay demand Z_r ----------------------
% Z_r = sum over INF stations of V(i,r) * mean service time(i,r).
ST = 1 ./ sn.rates;
ST(~isfinite(ST)) = 0;
Z = zeros(1, K);
for ist = find(isINF(:))'
    for r = 1:K
        Z(r) = Z(r) + V(ist, r) * ST(ist, r);
    end
end
D = zeros(M, K);                          % per-class demand at BCMP queues
for ist = find(isQ(:))'
    for r = 1:K
        D(ist, r) = V(ist, r) * ST(ist, r);
    end
end

% OI-station class visit ratios feed the v-weighted balanced-fairness balance
% (pfqn_ncoi / oi_phi); general (non-unit) visits are supported.
oiList = find(isOI(:))';
oivis = cell(1, numel(oiList));
for m = 1:numel(oiList)
    oivis{m} = V(oiList(m), :);
end

% ---- OI rank-rate handles on a per-class count vector ----------------------
% svcRateFun(c) is permutation-invariant, so it is a function of the count vector
% n; evaluate it on a canonical microstate holding n_r copies of class r. see
% _kb/06-solver-catalog.md (NC section, OI analyzer)
rates = cell(1, numel(oiList));
for m = 1:numel(oiList)
    fun = svc{oiList(m)};
    rates{m} = @(n) fun(oi_microstate(n));
end

% ---- population lattice ----------------------------------------------------
[shp, stride, total] = oi_lattice(N);

% ---- core normalizing-constant table (OI stations + aggregated delay) ------
% One call: the balanced-fairness convolution is a lattice convolution, so its
% internal table already holds G(n) for every 0 <= n <= N on the same
% column-major stride used here. Re-calling it per population would cost a
% needless factor total = prod_r (N_r+1).
[~, ~, Gfull] = pfqn_ncoi(Z, N, rates, oivis);

% ---- fold the BCMP queueing stations by lattice convolution ----------------
qList = find(isQ(:))';
for ist = qList
    c = sn.nservers(ist);
    Wq = oi_ld_table(D(ist, :), c, shp, total);
    Gfull = oi_conv(Gfull, Wq, shp, stride, total);
end

G = Gfull(total);
lG = log(G);

% ---- per-class throughput X_r = G(N - e_r)/G(N) ----------------------------
X = zeros(1, K);
for r = 1:K
    if N(r) > 0
        er = zeros(1, K); er(r) = 1;
        X(r) = Gfull(1 + sum((N - er) .* stride)) / G;
    end
end

% ---- per-station per-class mean queue length via the FNC identity ----------
Q = zeros(M, K);
for m = 1:numel(oiList)                   % OI stations
    ist = oiList(m);
    Phi = oi_phi(rates{m}, N, oivis{m});
    for r = 1:K
        if N(r) > 0
            [~, Psir] = pfqn_oi_fnc(Phi, N, @(n) n(r));
            Q(ist, r) = oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1;
        end
    end
end
for ist = qList                           % BCMP queueing stations
    Wq = oi_ld_table(D(ist, :), sn.nservers(ist), shp, total);
    for r = 1:K
        if N(r) > 0
            [~, Psir] = pfqn_oi_fnc(Wq, N, @(n) n(r));
            Q(ist, r) = oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1;
        end
    end
end
for ist = find(isINF(:))'                 % delay: Little's law
    for r = 1:K
        Q(ist, r) = X(r) * V(ist, r) * ST(ist, r);
    end
end

% ---- throughput, utilization, response time per station --------------------
T = zeros(M, K);
U = zeros(M, K);
R = zeros(M, K);
for ist = 1:M
    for r = 1:K
        T(ist, r) = X(r) * V(ist, r);
    end
end
for ist = find(isINF(:))'
    U(ist, :) = Q(ist, :);                % INF utilization convention
end
for ist = qList                           % BCMP queue: offered-load per server
    c = sn.nservers(ist);
    if ~isfinite(c) || c <= 0, c = 1; end
    for r = 1:K
        U(ist, r) = X(r) * D(ist, r) / c;
    end
end
for m = 1:numel(oiList)
    % In-service utilization U_r = E[sir_r]/c via the functional-server identity
    % E[f(n)] = G^{+}/G - 1 (pfqn_oi_fnc, pfqn_oi_insvc); see
    % _kb/06-solver-catalog.md (NC section, OI analyzer)
    ist = oiList(m);
    S = sn.nservers(ist);
    if ~isfinite(S) || S <= 0, S = 1; end
    Phi = oi_phi(rates{m}, N, oivis{m});
    gins = pfqn_oi_insvc(rates{m}, N);
    for r = 1:K
        if N(r) > 0
            fr = @(n) gins(1 + sum(n .* stride), r);
            [~, Psir] = pfqn_oi_fnc(Phi, N, fr);
            U(ist, r) = (oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1) / S;
        end
    end
end
for ist = 1:M
    for r = 1:K
        if T(ist, r) > 0
            R(ist, r) = Q(ist, r) / T(ist, r);
        end
    end
end

C = zeros(1, K);                          % per-class system response time
for r = 1:K
    if X(r) > 0
        C(r) = N(r) / X(r);
    end
end

runtime = toc(Tstart);
end

% ==========================================================================
function [shp, stride, total] = oi_lattice(N)
% Column-major lattice descriptor for populations 0 <= n <= N.
N = round(N(:)'); R = numel(N); shp = N + 1;
stride = ones(1, R);
for d = 2:R, stride(d) = stride(d-1) * shp(d-1); end
total = prod(shp);
end

% ==========================================================================
function c = oi_microstate(n)
% Canonical ordered microstate holding n_r copies of class r (n a count
% vector). For an order-independent station the service rate is invariant to
% the ordering, so this representative suffices to evaluate svcRateFun(c).
c = repelem(1:numel(n), round(n));
end

% ==========================================================================
function n = oi_sub(i, shp)
% Decode linear index i (1-based) to the subscript vector n (0-based counts).
R = numel(shp); n = zeros(1, R); li = i - 1;
for d = 1:R, n(d) = mod(li, shp(d)); li = floor(li / shp(d)); end
end

% ==========================================================================
function Phi = oi_phi(oirate, N, vis)
% Forward v-weighted balanced-fairness fill of the OI balance function:
% Phi(0)=1, Phi(n) = (1/mu(n)) sum_{r: n_r>0} v_r Phi(n - e_r).
[shp, stride, total] = oi_lattice(N);
R = numel(shp);
if nargin < 3 || isempty(vis), vis = ones(1, R); end
Phiv = zeros(total, 1);
for i = 1:total
    n = oi_sub(i, shp);
    if sum(n) == 0, Phiv(i) = 1; continue, end
    s = 0;
    for r = 1:R, if n(r) > 0, s = s + vis(r) * Phiv(i - stride(r)); end, end
    Phiv(i) = s / oirate(n);
end
if R == 1, Phi = Phiv; else, Phi = reshape(Phiv, shp); end
end

% ==========================================================================
function W = oi_ld_table(Dq, c, shp, total)
% BCMP load-dependent weight table over the lattice:
%   W(n) = (sum n)!/prod(n_r!) * prod_r D_r^{n_r} / prod_{k=1}^{sum n} beta(k),
% beta(k) = min(k,c) for a c-server queue (c=1 -> single server, beta==1).
% Column-major flat vector. W(0)=1.
Dq = Dq(:)'; R = numel(shp); W = zeros(total, 1);
if ~isfinite(c) || c <= 0, c = 1; end
for i = 1:total
    n = oi_sub(i, shp);
    tot = sum(n);
    logf = gammaln(tot + 1); ok = true;
    for r = 1:R
        if n(r) > 0
            if Dq(r) <= 0, ok = false; break, end
            logf = logf + n(r) * log(Dq(r)) - gammaln(n(r) + 1);
        end
    end
    if ~ok, continue, end
    for k = 1:tot
        logf = logf - log(min(k, c));
    end
    W(i) = exp(logf);
end
end

% ==========================================================================
function Cv = oi_conv(Av, Bv, shp, stride, total)
% Lattice convolution Cv(m) = sum_{0<=a<=m} Av(a) Bv(m-a) over 0..N.
R = numel(shp);
subs = zeros(total, R);
for i = 1:total, subs(i, :) = oi_sub(i, shp); end
Cv = zeros(total, 1);
for i = 1:total
    m = subs(i, :); acc = 0;
    for j = 1:i
        a = subs(j, :);
        if all(a <= m)
            acc = acc + Av(j) * Bv(1 + sum((m - a) .* stride));
        end
    end
    Cv(i) = acc;
end
end

% ==========================================================================
function val = oi_fnc_mean(Psi, Gfull, shp, stride, total)
% G^{+} = sum_{0<=b<=N} Psi(b) G(N-b): the FNC of the target station convolved
% against the full-network normalizing-constant table, evaluated at n = N.
N = shp - 1; Psiv = Psi(:); val = 0;
for i = 1:total
    if Psiv(i) == 0, continue, end
    b = oi_sub(i, shp);
    val = val + Psiv(i) * Gfull(1 + sum((N - b) .* stride));
end
end
