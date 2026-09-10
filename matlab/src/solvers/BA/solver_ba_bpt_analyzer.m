function [Q,U,R,T,C,X,lG,runtime] = solver_ba_bpt_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_BA_BPT_ANALYZER(SN, OPTIONS)
%
% Achievable-region LOWER bound on the mean response times of a multiclass
% open Markovian network, valid for EVERY non-idling scheduling policy at
% every station. The polyhedron is the first-order linear-programming
% relaxation of the achievable region (NPFQN_BND_BPT); this analyzer maps the
% LINE model onto it and reads the bound back per station and class.
%
% CLASS SPACE. The reference's "class" is a buffer: one exponential service
% rate, one Markovian routing law. LINE's (station, job class) pair is exactly
% that, so a pair carrying traffic becomes one LP class, the Source is absorbed
% into the external arrival vector, and class switching needs no special
% treatment because sn.rt already carries it.
%
% BOUND CONVENTION. R(i,r) is obtained by minimizing x over the polyhedron
% with the objective set to the unit vector of that pair, so each entry is a
% valid lower bound on its own. Q follows by Little's law from the bounded
% R and the EXACT throughput T (an open network's per-class rates are fixed by
% the traffic equations, not by the policy), and so does the per-class system
% response time C. U is exact for the same reason.
%
% TIGHTNESS. The relaxation is exact on M/M/1 and tight on the externally fed
% classes, but it is weak on a class whose arrivals are all internal: the
% only term coupling x_r to the second-moment block carries the factor
% lambda0_r, so an internally fed class can fall back to its own mean service
% time. That is a property of the first-order relaxation, not of this port;
% the reference's remedy is the higher-order (nonlinear/semidefinite)
% characterizations of its Section 5.
%
% Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994). Optimization
% of multiclass queueing networks: polyhedral and nonlinear characterizations
% of achievable performance. Annals of Applied Probability 4(1), 43-75.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;

Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K); lG = NaN;

% ----- model gates -----
% BA_OPEN_REFUSAL is the predicate SolverBA.supportsModelMethod reports, so the
% run raises the sentence the report gave; the checks below are its guards.
reason = ba_open_refusal(sn, 'bpt.lower');
if ~isempty(reason)
    line_error(mfilename, '%s', reason);
end
if any(isfinite(sn.njobs))
    line_error(mfilename, ...
        'Method ''bpt.lower'' supports fully open networks only (no closed classes).');
end
isSource = false(M,1);
for i = 1:M
    isSource(i) = (sn.nodetype(sn.stationToNode(i)) == NodeType.Source);
end
if ~any(isSource)
    line_error(mfilename, 'Method ''bpt.lower'' requires an open network with a Source station.');
end
qstat = find(~isSource);
if any(sn.sched(qstat) == SchedStrategy.INF)
    line_error(mfilename, ...
        'Method ''bpt.lower'' does not support delay (infinite-server) stations: the achievable region is derived for one server per station.');
end
if any(sn.nservers(qstat) > 1)
    line_error(mfilename, ...
        'Method ''bpt.lower'' does not support multi-server stations.');
end

% ----- station-space routing, with the Source absorbed into lambda0 -----
rtst = sn_rt_stations(sn);

% Pair (station i, class r) -> flat index (i-1)*K + r in the station space.
nq = numel(qstat);
pairStation = zeros(nq*K,1);
pairClass = zeros(nq*K,1);
pairFlat = zeros(nq*K,1);
np = 0;
for a = 1:nq
    i = qstat(a);
    for r = 1:K
        np = np + 1;
        pairStation(np) = i;
        pairClass(np) = r;
        pairFlat(np) = (i-1)*K + r;
    end
end

% External injection: whatever the Source stations emit, routed one hop.
lambda0 = zeros(np,1);
srcList = find(isSource);
for si = 1:numel(srcList)
    s = srcList(si);
    for r0 = 1:K
        arr = sn.rates(s,r0);
        if ~isfinite(arr) || arr <= 0
            continue
        end
        srow = (s-1)*K + r0;
        for p = 1:np
            lambda0(p) = lambda0(p) + arr * rtst(srow, pairFlat(p));
        end
    end
end

% Pair-to-pair routing. Flow to the Sink or back to a Source is the exit
% probability, i.e. the row deficit, and needs no column.
P = zeros(np,np);
for p = 1:np
    for q = 1:np
        P(p,q) = rtst(pairFlat(p), pairFlat(q));
    end
end

% ----- restrict to the pairs that actually carry traffic -----
lamAll = (eye(np) - P') \ lambda0;
keep = find(lamAll > 1e-12 * max(1, max(lamAll)));
if isempty(keep)
    line_error(mfilename, 'The model carries no open traffic.');
end
lambda0 = lambda0(keep);
P = P(keep,keep);
pairStation = pairStation(keep);
pairClass = pairClass(keep);
np = numel(keep);

mu = zeros(np,1);
for p = 1:np
    mu(p) = sn.rates(pairStation(p), pairClass(p));
    if ~isfinite(mu(p)) || mu(p) <= 0
        line_error(mfilename, ...
            'Station %d has no service rate for class %d but carries its traffic.', ...
            pairStation(p), pairClass(p));
    end
    if sn.procid(pairStation(p), pairClass(p)) ~= ProcessType.EXP
        line_error(mfilename, ...
            'Method ''bpt.lower'' requires exponential service: station %d class %d is %s.', ...
            pairStation(p), pairClass(p), ProcessType.toText(sn.procid(pairStation(p), pairClass(p))));
    end
end

% Renumber the stations the LP sees so its station index space is dense.
[~, ~, stationOf] = unique(pairStation);

% ----- one LP per pair, objective = that pair's unit vector -----
for p = 1:np
    e = zeros(np,1); e(p) = 1;
    [zlb, ~, info] = npfqn_bnd_bpt(lambda0, mu, P, stationOf, e);
    i = pairStation(p); r = pairClass(p);
    R(i,r) = zlb;
    T(i,r) = info.lambda(p);
    U(i,r) = info.rho(p);
end

% ----- exact open-network quantities -----
for si = 1:numel(srcList)
    s = srcList(si);
    for r = 1:K
        arr = sn.rates(s,r);
        if isfinite(arr) && arr > 0
            T(s,r) = T(s,r) + arr;
            X(1,r) = X(1,r) + arr;
        end
    end
end
Q = T .* R;
for r = 1:K
    if X(1,r) > 0
        C(1,r) = sum(Q(:,r)) / X(1,r);
    end
end
runtime = toc(T0);
end
