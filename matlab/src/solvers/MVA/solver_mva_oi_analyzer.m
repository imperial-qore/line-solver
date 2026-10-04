function [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_mva_oi_analyzer(sn, options)
% SOLVER_MVA_OI_ANALYZER Exact mean-value MVA for order-independent networks
%
% An order-independent (OI) station is a class-dependent load-dependent server
% whose total service rate mu(n) is a permutation-invariant function of the
% per-class count vector n. A closed network of infinite-server (delay) and
% load-independent (single-server, product-form) stations plus ANY number of OI
% stations is product-form. This analyzer aggregates the delay stations into a
% single think-time vector Z, collects the load-independent (LI) queue demands
% and the OI-station rate handles, and calls PFQN_MVAOI, the mean-value
% Conditional-MVA (CMVA) that carries one rate-shift vector per OI station and
% returns exact per-class throughput and queue-lengths WITHOUT any normalizing
% constant or joint marginal. The marginal-distribution counterpart is
% PFQN_MVAOI_MARG.
%
% Reference:
%   Reiser, Lavenberg (1980). Mean-Value Analysis of Closed Multichain Queuing
%   Networks. JACM 27(2). Load-dependent extension: Bruell, Balbo, Afshari
%   (1984). OI stations / CMVA: Casale (2009); Casale, Comte, Dorsman (2026).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tstart = tic;

oi_list = find_oi_stations(sn);
if isempty(oi_list)
    error('solver_mva_oi_analyzer:NoOIStation', ...
        'OI solver requires at least one order-independent station');
end

% ---- reject class switching (OI rank rates are per raw class) --------------
% The recursion is driven by the per-class population vector sn.njobs, which
% class switching makes meaningless: a class that only ever appears mid-chain
% carries njobs = 0, so the OI station would be analyzed as if empty. Refuse it
% the way SOLVER_NC_OI_ANALYZER does rather than return that silently.
for c = 1:sn.nchains
    if numel(sn.inchain{c}) > 1
        line_error(mfilename, 'solver_mva_oi requires one class per chain (no class switching).');
    end
end

M = sn.nstations;
R = sn.nclasses;
N = round(sn.njobs(:)');

isOI = false(1, M);
isOI(oi_list) = true;
isDelay = (sn.sched(:)' == SchedStrategy.INF);

% Per-class service demand D_ir = V_ir/rate_ir at non-OI stations.
V = zeros(M, R);
for c = 1:numel(sn.visits)
    V = V + sn.visits{c};
end
D = zeros(M, R);
for i = 1:M
    for r = 1:R
        if isfinite(sn.rates(i,r)) && sn.rates(i,r) > 0
            D(i,r) = V(i,r) / sn.rates(i,r);
        end
    end
end

% Aggregate the delay stations into Z; collect the LI queue demands and OI-station
% rate handles. A multiserver BCMP queue (c>1) cannot take the single-server LI
% path, so it is promoted to an OI station via the c-server BCMP weight; see
% _kb/06-solver-catalog.md (NC section, OI analyzer) for the W(n)/mu(n) identity.
Z = zeros(1, R);
li_list = [];
ms_list = [];
for i = 1:M
    if isOI(i)
        continue
    elseif isDelay(i)
        Z = Z + D(i,:);
    elseif isfinite(sn.nservers(i)) && sn.nservers(i) > 1
        ms_list(end+1) = i; %#ok<AGROW>
    else
        li_list(end+1) = i; %#ok<AGROW>
    end
end
Dli = D(li_list, :);

muCell = cell(1, numel(oi_list) + numel(ms_list));
% Per-muCell-station visit vectors. Genuine OI stations carry their class visit
% ratios V(oi,:) (the rate handle has no visits); ms-promoted stations pass
% ones, as their visits are already folded into D by ms_oi_rate.
oivis = cell(1, numel(oi_list) + numel(ms_list));
for o = 1:numel(oi_list)
    node_oi = sn.stationToNode(oi_list(o));
    svcRateFun = sn.nodeparam{node_oi}.svcRateFun;
    muCell{o} = @(n) oi_rate(svcRateFun, n, R);
    oivis{o} = V(oi_list(o), :);
end
for j = 1:numel(ms_list)
    muCell{numel(oi_list) + j} = ms_oi_rate(D(ms_list(j), :), sn.nservers(ms_list(j)));
    oivis{numel(oi_list) + j} = ones(1, R);
end

[X, Qoi, Qli, ~, Soi] = pfqn_mvaoi(Z, N, muCell, Dli, oivis, options);

% Assemble per-station mean queue-lengths.
QN = zeros(M, R);
for o = 1:numel(oi_list)
    QN(oi_list(o), :) = Qoi(o, :);
end
for j = 1:numel(ms_list)
    QN(ms_list(j), :) = Qoi(numel(oi_list) + j, :);
end
for j = 1:numel(li_list)
    QN(li_list(j), :) = Qli(j, :);
end
for i = 1:M
    if isDelay(i)
        QN(i, :) = X .* D(i, :);     % IS station: exact product-form share
    end
end
XN = X;

% Row of Soi/Qoi holding each OI station (muCell order: oi_list, then ms_list).
oiRow = zeros(1, M);
for o = 1:numel(oi_list)
    oiRow(oi_list(o)) = o;
end

% Assemble outputs at the full population.
TN = zeros(M, R);
RN = zeros(M, R);
UN = zeros(M, R);
CN = zeros(M, R);
for i = 1:M
    for r = 1:R
        TN(i,r) = XN(r) * V(i,r);
        if XN(r) > 0
            RN(i,r) = QN(i,r) / XN(r);
        end
        if isOI(i)
            % In-service utilization U_r = E[sir_r]/nservers; see
            % _kb/06-solver-catalog.md (NC section, OI analyzer). Soi holds E[sir_r].
            sv = sn.nservers(i);
            if ~isfinite(sv) || sv <= 0
                sv = 1;
            end
            UN(i,r) = Soi(oiRow(i), r) / sv;
        elseif isDelay(i)
            UN(i,r) = QN(i,r);
        else
            sv = sn.nservers(i);
            if ~isfinite(sv) || sv <= 0
                sv = 1;
            end
            UN(i,r) = XN(r) * D(i,r) / sv;
        end
        CN(i,r) = RN(i,r);
    end
end

lG = 0;
runtime = toc(tstart);
iter = sum(N);
actualmethod = 'oi';
end

% =========================================================================
% Helper functions
% =========================================================================

function oi_list = find_oi_stations(sn)
% Station indices of all OI stations: PAS/OI scheduling with an all-zero swap
% graph and a service-rate function (mirrors nc_is_oi_model detection).
oi_list = [];
for ist = 1:sn.nstations
    if sn.sched(ist) ~= SchedStrategy.PAS && sn.sched(ist) ~= SchedStrategy.OI
        continue
    end
    ind = sn.stationToNode(ist);
    if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind})
        continue
    end
    np = sn.nodeparam{ind};
    if ~isfield(np, 'swapGraph') || ~isfield(np, 'svcRateFun') || isempty(np.svcRateFun)
        continue
    end
    sg = np.swapGraph;
    if isempty(sg) || any(sg(:) ~= 0)
        continue
    end
    oi_list(end+1) = ist; %#ok<AGROW>
end
end

function h = ms_oi_rate(Dq, c)
% OI rate function reproducing the c-server BCMP station with per-class demands
% Dq: mu(n) = (min(|n|,c)/|n|) * sum_{r: n_r>0} n_r/Dq_r. For c = 1 this is the
% familiar total completion rate of a multiclass single-server queue, and for a
% single class it reduces to min(n,c)/Dq (M/M/c).
Dq = Dq(:)';
if ~isfinite(c) || c <= 0
    c = 1;
end
h = @(n) ms_oi_rate_eval(n, Dq, c);
end

function rate = ms_oi_rate_eval(n, Dq, c)
tot = sum(n);
if tot == 0
    rate = 0;
    return
end
acc = 0;
for r = 1:numel(n)
    if n(r) > 0 && Dq(r) > 0
        acc = acc + n(r) / Dq(r);
    end
end
rate = (min(tot, c) / tot) * acc;
end

function rate = oi_rate(svcRateFun, n, R)
% OI total service rate at count vector n via the 1-based canonical microstate.
if sum(n) == 0
    rate = 0;
    return
end
micro = repelem(1:R, n);
rate = svcRateFun(micro);
end
