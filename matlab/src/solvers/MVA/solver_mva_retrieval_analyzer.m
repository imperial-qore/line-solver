function [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,delayedprob,hitproblist,itemprob,latency,runtime,method] = solver_mva_retrieval_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,HITPROB,MISSPROB,DELAYEDPROB,HITPROBLIST,LATENCY,RUNTIME,METHOD] = SOLVER_MVA_RETRIEVAL_ANALYZER(SN, OPTIONS)
%
% Fixed-point approximation of a delayed-hit (retrieval-system) cache via the FPI
% algorithms: retrieval_fpi (hit / miss / delayed-hit ratios) and
% retrieval_fpi_latency (expected latency Z, eq:latency tot).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
CN = [];
XN = zeros(1, sn.nclasses);
QN = zeros(sn.nstations, sn.nclasses);
UN = zeros(sn.nstations, sn.nclasses);
RN = zeros(sn.nstations, sn.nclasses);
TN = zeros(sn.nstations, sn.nclasses);
lG = NaN;
method = 'fpi';

line_debug('MVA retrieval analyzer starting: nclasses=%d', sn.nclasses);

[m, lambda, gamma, eta, alpha, T, R, station_type] = cache_retrieval_inputs(sn);
n = numel(lambda);

% --- FPI hit / miss / delayed-hit ratios ---
[pmiss, phit, pdh] = retrieval_fpi(m, lambda, eta, gamma);
pi0  = pmiss(:).';      % per-item miss
pih  = sum(phit, 1);    % per-item hit
phid = sum(pdh, 1);     % per-item delayed

% --- FPI expected latency Z (eq:latency tot) ---
Z = retrieval_fpi_latency(m, lambda, gamma, alpha, T, R, station_type);

% --- read class and per-class aggregates (single read class) ---
ci = find(sn.nodetype == NodeType.Cache);
ch = sn.nodeparam{ci};
rk = keys(ch.retrievalSystemQueueIndices);
jobinClass = double(rk{1}) + 1;
w = lambda(:) / sum(lambda);
hitAgg = sum(w .* pih(:));
missAgg = sum(w .* pi0(:));
delayedAgg = sum(w .* phid(:));

% --- source throughput ---
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
sourceRate = sn.rates(source_ist, :);
sourceRate(isnan(sourceRate)) = 0;
TN(source_ist, :) = sourceRate;

% --- cache results (hitprob is the TRUE hit fraction; delayedprob the
% delayed-hit fraction; true hit + delayed + miss = 1) ---
hitprob     = NaN(1, sn.nclasses);
missprob    = NaN(1, sn.nclasses);
delayedprob = NaN(1, sn.nclasses);
latency     = NaN(1, sn.nclasses);
hitprob(jobinClass)     = hitAgg;
missprob(jobinClass)    = missAgg;
delayedprob(jobinClass) = delayedAgg;
latency(jobinClass)     = Z;

% per-list (per-level) hit fractions for the read class: phit is (h x n);
% access-weighted over items gives the per-list hit probability (rows sum to
% hitAgg over lists).
hitproblist = NaN(sn.nclasses, numel(m));
hitproblist(jobinClass, :) = (phit * w(:)).';

% per-item occupancy [nitems x (lists+1)]: column 1 = miss (item not cached),
% columns 2..end = probability the item resides in each list.
itemprob = [pi0(:), phit.'];

hc = ch.hitclass(jobinClass);
mc = ch.missclass(jobinClass);
if hc > 0, XN(hc) = sourceRate(jobinClass) * (hitAgg + delayedAgg); end
if mc > 0, XN(mc) = sourceRate(jobinClass) * missAgg; end

% --- retrieval-station mean occupancy (QLen) and throughput ---
queueNodes = double(ch.retrievalSystemQueueIndices(rk{1}));
S = numel(queueNodes);
psIdx = find(station_type == "PS" | station_type == "SIRO" | station_type == "FCFS" | station_type == "LCFSPR");   % SIRO/FCFS/LCFSPR as PS
for s = 1:S
    sst = sn.nodeToStation(queueNodes(s));
    if station_type(s) == "PS" || station_type(s) == "SIRO" || station_type(s) == "FCFS" || station_type(s) == "LCFSPR"
        prow = 1 + find(psIdx == s);
    else
        prow = 1;
    end
    phi_s = 0;
    if prow <= size(pdh, 1)
        phi_s = sum(pdh(prow, :));
    end
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
        RN(sst, jobinClass) = phi_s / tput_s;
    end
end

runtime = toc(T0);
end
