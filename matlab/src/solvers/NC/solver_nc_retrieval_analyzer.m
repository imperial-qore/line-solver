function [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,delayedprob,hitproblist,itemprob,latency,runtime,method] = solver_nc_retrieval_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,HITPROB,MISSPROB,DELAYEDPROB,HITPROBLIST,LATENCY,RUNTIME,METHOD] = SOLVER_NC_RETRIEVAL_ANALYZER(SN, OPTIONS)
%
% Exact analysis of a delayed-hit (retrieval-system) cache via the product-form
% recurrence algorithms: retrieval_nc (normalizing constant) and retrieval_metrics
% (hit / miss / delayed-hit ratios). Latency is left to SolverMVA
% (retrieval_fpi_latency) and returned as NaN here.
%
% Those recurrences are exponential in the number of items, so OPTIONS.METHOD =
% 'rayint' selects instead the ray (WKB) approximation of retrieval_rayint, which
% is polynomial. It applies only when every fetch station is infinite-server,
% where the delayed-hit constant factorizes exactly as prod_k D_k times the plain
% cache constant with access factors gamma_{k,j}/D_k; elsewhere it warns and falls
% back to the exact path. Accuracy is that of retrieval_rayint, i.e. roughly
% 0.14*(1/min_j m_j + 1/(n - sum_j m_j)) in relative terms.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
CN = [];
XN = zeros(1, sn.nclasses);
QN = zeros(sn.nstations, sn.nclasses);
UN = zeros(sn.nstations, sn.nclasses);
RN = zeros(sn.nstations, sn.nclasses);
TN = zeros(sn.nstations, sn.nclasses);
method = 'exact';

line_debug('NC retrieval analyzer starting: nclasses=%d', sn.nclasses);

[m, lambda, gamma, eta, alpha, T, R, station_type] = cache_retrieval_inputs(sn); %#ok<ASGLU>
n = numel(lambda);
r = sum(station_type == "PS");

useray = isfield(options, 'method') && any(strcmpi(options.method, {'rayint','ray'}));
if useray
    % The ray expansion needs the delayed-hit constant to factorize. It does,
    % EXACTLY, when every fetch station is infinite-server: dividing the
    % retrieval_nc recurrence by prod_k D_k with D_k = 1 + lambda_k eta_{0,k}
    % collapses it onto cache_erec with theta_{k,j} = gamma_{k,j}/D_k, so the
    % delayed-hit cache IS a plain cache with fetch-inflated access factors.
    % A queueing (PS) fetch station breaks this: the (v_s+1) multiplicity ties
    % E(0,m) to the whole moment tower E(1_s,m), E(2_s,m), ..., and replacing it
    % by the retrieval_fpi mean field overestimates E by 13%/140%/830% at
    % n=6/8/10 (measured), growing with n. There is no validated expansion
    % there, so refuse rather than return a confident wrong number.
    reason = '';
    if size(eta, 2) > 1 && any(any(eta(:, 2:end) ~= 0))
        reason = 'the retrieval system has a queueing (non infinite-server) fetch station';
    elseif sum(m) >= n
        reason = 'the cache is full (sum(m) >= n), where the saddle point escapes to infinity';
    end
    if ~isempty(reason)
        line_warning(mfilename, ['method ''rayint'' does not apply because %s; ' ...
            'falling back to the exact recurrences.\n'], reason);
        useray = false;
    end
end

if useray
    % --- ray (WKB) approximation, infinite-server fetch ---
    D     = 1 + lambda(:).*eta(:, 1);
    theta = gamma ./ D;
    [~, lGcache, rayout] = retrieval_rayint(theta, m);
    lG = sum(log(D)) + lGcache;

    % Same saddle as the constant, so the ratios are consistent with lG:
    % pi_{i,j} = theta_{i,j} xi_j / (1 + sum_l theta_{i,l} xi_l), and the
    % out-of-cache mass 1 - sum_j pi_{i,j} splits between a true miss (weight 1)
    % and an outstanding fetch (weight lambda_i eta_{0,i}) in proportion 1:D_i-1.
    xi   = rayout.xi(:).';
    txi  = theta .* xi;
    phit = (txi ./ (1 + sum(txi, 2))).';        % h x n
    pmiss = ((1 - sum(phit, 1)).') ./ D;        % pi_{i,0} = (1 - pihit_i)/D_i
    pmiss = pmiss(:).';
    pdh   = (lambda(:).*eta(:, 1)).' .* pmiss;  % phi_{0,i}, single IS row
    method = 'rayint';
else
    % --- exact normalizing constant E(m) = retrieval_nc(0,m,...) ---
    E = retrieval_nc(zeros(1, r), m, lambda, eta, gamma);
    lG = log(E);

    % --- exact hit / miss / delayed-hit ratios ---
    [pmiss, phit, pdh] = retrieval_metrics(m, lambda, eta, gamma);
end
pi0  = pmiss(:).';      % per-item miss   pi_{i,0}
pih  = sum(phit, 1);    % per-item hit    sum_j pi_{i,j}
phid = sum(pdh, 1);     % per-item delayed sum_s phi_{s,i}

% --- read class and per-class aggregates (single read class) ---
ci = find(sn.nodetype == NodeType.Cache);
ch = sn.nodeparam{ci};
rk = keys(ch.retrievalSystemQueueIndices);
jobinClass = double(rk(1)) + 1;
w = lambda(:) / sum(lambda);            % access-weighted item mixture
hitAgg = sum(w .* pih(:));
missAgg = sum(w .* pi0(:));
delayedAgg = sum(w .* phid(:));

% --- source throughput ---
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
sourceRate = sn.rates(source_ist, :);
sourceRate(isnan(sourceRate)) = 0;
TN(source_ist, :) = sourceRate;

% --- cache results (per-class vectors; read-class entry set, rest NaN).
% hitprob is the TRUE hit fraction; delayedprob the delayed-hit fraction;
% true hit + delayed + miss = 1 ---
hitprob     = NaN(1, sn.nclasses);
missprob    = NaN(1, sn.nclasses);
delayedprob = NaN(1, sn.nclasses);
latency     = NaN(1, sn.nclasses);
hitprob(jobinClass)     = hitAgg;
missprob(jobinClass)    = missAgg;
delayedprob(jobinClass) = delayedAgg;  % delayed implicit: 1 - hit - miss = sum_s phi

% per-list (per-level) hit fractions for the read class: phit is (h x n);
% access-weighted over items gives the per-list hit probability.
hitproblist = NaN(sn.nclasses, numel(m));
hitproblist(jobinClass, :) = (phit * w(:)).';

% per-item occupancy [nitems x (lists+1)]: column 1 = miss (item not cached),
% columns 2..end = probability the item resides in each list.
itemprob = [pi0(:), phit.'];

% --- throughputs: misses -> missclass, hits+delayed -> hitclass ---
hc = ch.hitclass(jobinClass);
mc = ch.missclass(jobinClass);
if hc > 0, XN(hc) = sourceRate(jobinClass) * (hitAgg + delayedAgg); end
if mc > 0, XN(mc) = sourceRate(jobinClass) * missAgg; end

% --- retrieval-station mean occupancy (QLen) and throughput from phi/visits ---
% phi_{s,i} (pdh) is the mean number of item i being retrieved at station s;
% summing over items gives the station occupancy.
queueNodes = double(ch.retrievalSystemQueueIndices{rk(1)});
S = numel(queueNodes);
isIdx = find(station_type == "IS"); %#ok<NASGU>
psIdx = find(station_type == "PS" | station_type == "SIRO" | station_type == "FCFS" | station_type == "LCFSPR");   % SIRO/FCFS/LCFSPR as PS
for s = 1:S
    sst = sn.nodeToStation(queueNodes(s));
    if station_type(s) == "PS" || station_type(s) == "SIRO" || station_type(s) == "FCFS" || station_type(s) == "LCFSPR"
        prow = 1 + find(psIdx == s);    % pdh row for this PS/SIRO/FCFS/LCFSPR station
    else
        prow = 1;                       % IS aggregate row
    end
    phi_s = 0;
    if prow <= size(pdh, 1)
        phi_s = sum(pdh(prow, :));      % station occupancy = utilization (paper phi_s)
    end
    % fetch throughput through station s = sum_i (item fetch rate)*visits_{s,i}
    tput_s = 0;
    for i = 1:n
        Ri = R(:, :, i);
        a = Ri(1, 2:S+1);
        Pmat = Ri(2:S+1, 2:S+1);
        vis = a / (eye(S) - Pmat);
        tput_s = tput_s + sourceRate(jobinClass) * (lambda(i)/sum(lambda)) * pi0(i) * vis(s);
    end
    QN(sst, jobinClass) = phi_s;
    UN(sst, jobinClass) = phi_s;
    TN(sst, jobinClass) = tput_s;
    if tput_s > 0
        RN(sst, jobinClass) = phi_s / tput_s;   % Little's law sojourn time
    end
end

runtime = toc(T0);
end
