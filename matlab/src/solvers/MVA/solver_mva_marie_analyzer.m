function [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_mva_marie_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,LG,RUNTIME,LASTITER,ACTUALMETHOD] = SOLVER_MVA_MARIE_ANALYZER(SN, OPTIONS)
%
% Marie's iterative aggregation-decomposition (Marie 1979/1980) for closed
% networks with FCFS non-exponential (Coxian) service, wired as SolverMVA
% method 'marie'. Infinite-server stations fold into per-chain think time;
% pfqn_marie is applied to the queueing stations, where only FCFS is service-
% sensitive (its SCV is used) while insensitive product-form disciplines
% (PS, LCFSPR) are forced to exponential. Chain results are then deaggregated
% to classes via sn_deaggregate_chain_results.
%
% Single chain: aggregate is the exact load-dependent product-form solve, so
% exponential service is exact. Multiple chains: aggregate is QD-AMVA with
% class-dependent scaling, exact for product form and approximate otherwise.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
lG = NaN;
actualmethod = 'marie';

M = sn.nstations;
[Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn); %#ok<ASGLU>
C = sn.nchains;

%% gating
% Closed models only.
if any(isinf(Nchain)) || any(sn.nodetype == NodeType.Source)
    line_error(mfilename, ['The ''marie'' method supports closed models only; this model has open ' ...
        'classes. Use another SolverMVA method (e.g. ''default'').']);
end

% Delay/infinite-server stations (fold into think time); the rest are queueing.
isDelay = isinf(sn.nservers(:)) | (sn.sched(:) == SchedStrategy.INF);

% Scheduling support: FCFS (service-sensitive), the insensitive product-form
% disciplines PS/LCFSPR (treated as exponential), and Delay. Reject others.
schedOK = isDelay | (sn.sched(:) == SchedStrategy.FCFS) | ...
    (sn.sched(:) == SchedStrategy.PS) | (sn.sched(:) == SchedStrategy.LCFSPR);
if ~all(schedOK)
    bad = find(~schedOK, 1);
    line_error(mfilename, sprintf(['The ''marie'' method supports FCFS, PS, LCFSPR and Delay ' ...
        'stations only; station %d has an unsupported scheduling strategy. Use another SolverMVA method.'], bad));
end

queueRows = find(~isDelay);
delayRows = find(isDelay);
Mq = numel(queueRows);

%% per-chain think time from the delay stations
Z = zeros(1,C);
for c = 1:C
    Z(c) = sum(Lchain(delayRows,c));
end

%% queueing-station demands and effective SCV
% Only FCFS is service-time sensitive; PS/LCFSPR are insensitive (product
% form), so their SCV is set to 1 (exponential-equivalent for Marie).
L = Lchain(queueRows,:);
SCV = SCVchain(queueRows,:);
for jj = 1:Mq
    if sn.sched(queueRows(jj)) ~= SchedStrategy.FCFS
        SCV(jj,:) = 1;
    end
end
SCV(~isfinite(SCV) | SCV <= 0) = 1;   % guard undefined SCV (e.g. zero demand)

% Multiserver isolation is supported for single-chain models only.
nservers = sn.nservers(queueRows);
nservers(~isfinite(nservers)) = 1;
if C > 1 && any(nservers > 1)
    line_error(mfilename, ['The ''marie'' method supports multiserver queueing stations for ' ...
        'single-chain models only; this model is multichain with a multiserver station.']);
end

%% Marie solve
if C == 1
    [Xm,Qm,Um,~,lastiter] = pfqn_marie(L, Nchain, Z, SCV, [], [], nservers);
else
    [Xm,Qm,Um,~,lastiter] = pfqn_marie(L, Nchain, Z, SCV);
end
Xchain = Xm(:)';                       % 1 x C chain throughput (reference station)

%% assemble full-station chain-level matrices
Qchain = zeros(M,C);
Uchain = zeros(M,C);
Rchain = zeros(M,C);
Qchain(queueRows,:) = Qm;
Uchain(queueRows,:) = Um;
Tchain = repmat(Xchain,M,1) .* Vchain; % per-station throughput = chain X * visits

% residence per station (Little's law, consistent with Marie's Q)
for ist = 1:M
    for c = 1:C
        if Tchain(ist,c) > 0
            Rchain(ist,c) = Qchain(ist,c) / Tchain(ist,c);
        end
    end
end

% delay stations: number in service Q = T*S, U = T*S, R = S
for jj = 1:numel(delayRows)
    ist = delayRows(jj);
    for c = 1:C
        Qchain(ist,c) = Tchain(ist,c) * STchain(ist,c);
        Uchain(ist,c) = Tchain(ist,c) * STchain(ist,c);
        Rchain(ist,c) = STchain(ist,c);
    end
end

[QN,UN,RN,TN,CN,XN] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, ...
    Qchain, Uchain, Rchain, Tchain, [], Xchain);

runtime = toc(Tstart);
end
