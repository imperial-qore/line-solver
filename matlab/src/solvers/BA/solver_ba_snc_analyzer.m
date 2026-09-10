function [Q,U,R,T,C,X,lG,runtime,env] = solver_ba_snc_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ENV] = SOLVER_BA_SNC_ANALYZER(SN, OPTIONS)
%
% Stochastic network calculus UPPER bound on the mean response times and queue
% lengths of a feed-forward open network, valid for EVERY work-conserving
% scheduling policy at every station. The api family in matlab/src/api/snc
% supplies the envelope algebra; this analyzer maps the LINE model onto it,
% propagates envelopes hop by hop, and reads the bound back per station and
% class.
%
% UNITS ARE JOBS, NOT WORK. The arrival envelope counts jobs and the service
% element is snc_srv_exp, the counting process of an Exp(mu) server. That is
% what lets a departure envelope from one station be the arrival envelope of
% the next: a service-time work unit differs from station to station, a job
% does not. On a single M/M/1 the resulting backlog bound decays as
% (lambda/mu)^n and the delay bound as exp(-(mu-lambda)*d), both exact rates.
%
% BOUND CONVENTION. R(i,r) is snc_mean_delay of the (arrival, service) envelope
% pair at that station, i.e. the integral of the delay tail bound, so each
% entry is a valid upper bound on its own. Q follows by Little's law from the
% bounded R and the EXACT throughput T (an open network's per-class rates are
% fixed by the traffic equations, not by the policy), and so does the per-class
% system response time C. U is exact for the same reason.
%
% WHAT IS ASSUMED, and refused when it does not hold:
%   - fully open network with a Source, no delay station, one server per station
%   - exponential service; the source may be any Markovian (D0,D1) process
%   - FEED-FORWARD at the station level: the station graph must be acyclic, so
%     that the cross traffic of a station is always determined upstream of it
%   - DETERMINISTIC ROUTING downstream of the Source: a probabilistic split of
%     an already-queued flow has no exponential-form envelope short of bounding
%     the split by the whole flow, which would report a false instability.
%     Splitting AT a Poisson Source is exact and is allowed
%   - classes sharing a station must have the SAME service rate, since the
%     blind-multiplexing leftover subtracts job counts from a job-count
%     capacity
%   - independence of the flow of interest and its cross traffic, and of the
%     stations; the dependent case needs a Hoelder split that the elementary
%     envelope algebra does not implement
%
% TIGHTNESS. This is a policy-robust bound, so it is loose on the mean: 2.4x
% the exact M/M/1 mean response time at rho = 0.1 and 10.4x at rho = 0.95. Its
% sharp object is the TAIL, whose decay rate it reproduces exactly; reach it
% through SolverBA.getDelayPerc / getBacklogPerc / getPercTable rather than
% through the mean columns when the quantile is what matters.
%
% Reference: M. Fidler and A. Rizk (2015). A Guide to the Stochastic Network
% Calculus. IEEE Communications Surveys and Tutorials 17(1), 92-105.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;
tol = 1e-10;

Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K); lG = NaN;
env = struct('arv', {cell(M,K)}, 'srv', {cell(M,K)}, 'carries', {false(M,K)});

% ----- model gates -----
% BA_OPEN_REFUSAL is the predicate SolverBA.supportsModelMethod reports, so the
% run raises the sentence the report gave (feed-forward graph, deterministic
% routing downstream of the Source, one rate per station, Poisson-only split);
% the checks below are its guards.
reason = ba_open_refusal(sn, 'snc.upper');
if ~isempty(reason)
    line_error(mfilename, '%s', reason);
end
if any(isfinite(sn.njobs))
    line_error(mfilename, ...
        'Method ''snc.upper'' supports fully open networks only (no closed classes).');
end
isSource = false(M,1);
for i = 1:M
    isSource(i) = (sn.nodetype(sn.stationToNode(i)) == NodeType.Source);
end
if ~any(isSource)
    line_error(mfilename, 'Method ''snc.upper'' requires an open network with a Source station.');
end
qstat = find(~isSource);
if any(sn.sched(qstat) == SchedStrategy.INF)
    line_error(mfilename, ...
        'Method ''snc.upper'' does not support delay (infinite-server) stations: the service envelope is that of a single busy server.');
end
if any(sn.nservers(qstat) > 1)
    line_error(mfilename, 'Method ''snc.upper'' does not support multi-server stations.');
end

% ----- station-space routing, with the Source absorbed into lambda0 -----
rtst = sn_rt_stations(sn);

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

% Source injection, kept per (source, class) so that the exogenous process of
% each stream is still identifiable once the pairs are known.
srcList = find(isSource);
inject = zeros(np, numel(srcList)*K);
lambda0 = zeros(np,1);
col = 0;
srcOfCol = zeros(1, numel(srcList)*K);
clsOfCol = zeros(1, numel(srcList)*K);
for si = 1:numel(srcList)
    s = srcList(si);
    for r0 = 1:K
        col = col + 1;
        srcOfCol(col) = s;
        clsOfCol(col) = r0;
        arr = sn.rates(s,r0);
        if ~isfinite(arr) || arr <= 0
            continue
        end
        srow = (s-1)*K + r0;
        for p = 1:np
            inject(p,col) = arr * rtst(srow, pairFlat(p));
            lambda0(p) = lambda0(p) + inject(p,col);
        end
    end
end

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
lam = lamAll(keep);
inject = inject(keep,:);
P = P(keep,keep);
pairStation = pairStation(keep);
pairClass = pairClass(keep);
np = numel(keep);

mu = zeros(np,1);
for p = 1:np
    i = pairStation(p); r = pairClass(p);
    mu(p) = sn.rates(i,r);
    if ~isfinite(mu(p)) || mu(p) <= 0
        line_error(mfilename, ...
            'Station %d has no service rate for class %d but carries its traffic.', i, r);
    end
    if sn.procid(i,r) ~= ProcessType.EXP
        line_error(mfilename, ...
            'Method ''snc.upper'' requires exponential service: station %d class %d is %s.', ...
            i, r, ProcessType.toText(sn.procid(i,r)));
    end
end

% ----- routing restrictions: no split downstream of the Source -----
for p = 1:np
    succ = find(P(p,:) > tol);
    if numel(succ) > 1
        line_error(mfilename, ...
            'Method ''snc.upper'' requires deterministic routing downstream of the Source: station %d class %d splits its flow over %d destinations.', ...
            pairStation(p), pairClass(p), numel(succ));
    end
    if numel(succ) == 1 && abs(P(p,succ) - 1) > 1e-8
        line_error(mfilename, ...
            'Method ''snc.upper'' requires deterministic routing downstream of the Source: station %d class %d routes onward with probability %g.', ...
            pairStation(p), pairClass(p), P(p,succ));
    end
end

% A Source that splits is exact only when its process is Poisson, since a
% Bernoulli thinning of a Poisson stream is again Poisson. Any other Markovian
% source must reach its station undivided.
for col = 1:size(inject,2)
    s = srcOfCol(col); r0 = clsOfCol(col);
    if s == 0 || sum(inject(:,col)) <= 0
        continue
    end
    dest = find(inject(:,col) > tol);
    if numel(dest) > 1 && sn.procid(s,r0) ~= ProcessType.EXP
        line_error(mfilename, ...
            'Method ''snc.upper'' can split only a Poisson Source: source %d class %d is %s and feeds %d stations.', ...
            s, r0, ProcessType.toText(sn.procid(s,r0)), numel(dest));
    end
end

% ----- one service rate per station, and a feed-forward station graph -----
stationsUsed = unique(pairStation);
for a = 1:numel(stationsUsed)
    i = stationsUsed(a);
    rates = mu(pairStation == i);
    if max(rates) - min(rates) > 1e-8 * max(1, max(rates))
        line_error(mfilename, ...
            'Method ''snc.upper'' requires the classes sharing a station to have equal service rates: station %d carries rates in [%g, %g].', ...
            i, min(rates), max(rates));
    end
end

ns = numel(stationsUsed);
stIdx = zeros(M,1);
stIdx(stationsUsed) = 1:ns;
adj = false(ns,ns);
for p = 1:np
    succ = find(P(p,:) > tol);
    if isempty(succ)
        continue
    end
    a = stIdx(pairStation(p));
    b = stIdx(pairStation(succ));
    adj(a,b) = true;
end
order = topoOrder(adj);
if isempty(order)
    line_error(mfilename, ...
        'Method ''snc.upper'' requires a feed-forward network: the station graph has a cycle, so a station''s cross traffic is not determined upstream of it.');
end

% ----- envelope propagation, station by station in feed-forward order -----
arvH = cell(np,1);
srvH = cell(np,1);
outH = cell(np,1);

for oi = 1:ns
    i = stationsUsed(order(oi));
    here = find(pairStation == i);

    % arrivals: exogenous injections plus the departure envelopes upstream
    for a = 1:numel(here)
        p = here(a);
        parts = {};
        for col = 1:size(inject,2)
            if inject(p,col) > tol
                s = srcOfCol(col); r0 = clsOfCol(col);
                if sn.procid(s,r0) == ProcessType.EXP
                    parts{end+1} = makePoisson(inject(p,col)); %#ok<AGROW>
                else
                    D = sn.proc{s}{r0};
                    parts{end+1} = makeMap(D{1}, D{2}); %#ok<AGROW>
                end
            end
        end
        pred = find(P(:,p)' > tol);
        for b = 1:numel(pred)
            parts{end+1} = outH{pred(b)}; %#ok<AGROW>
        end
        if isempty(parts)
            line_error(mfilename, ...
                'Station %d class %d carries traffic with no identifiable source.', i, pairClass(p));
        end
        arvH{p} = makeSum(parts);
    end

    % service: the busy-server counting process, minus the cross traffic that
    % blind multiplexing lets the other classes take from it
    for a = 1:numel(here)
        p = here(a);
        cross = cell(0,1);
        for b = 1:numel(here)
            if here(b) ~= p
                cross{end+1,1} = arvH{here(b)}; %#ok<AGROW>
            end
        end
        srvH{p} = makeLeftover(mu(p), cross);
        outH{p} = makeOutput(arvH{p}, srvH{p});
    end
end

% ----- read the bounds back -----
for p = 1:np
    i = pairStation(p); r = pairClass(p);
    R(i,r) = snc_mean_delay(arvH{p}, srvH{p});
    T(i,r) = lam(p);
    U(i,r) = lam(p) / mu(p);
    env.arv{i,r} = arvH{p};
    env.srv{i,r} = srvH{p};
    env.carries(i,r) = true;
end

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

% ----- envelope closures -----
% Each returns a handle of theta, so that the whole feed-forward composition is
% re-evaluated at whatever theta the Chernoff search asks for. Building the
% numbers eagerly at one theta would defeat the optimization.

function h = makePoisson(lambda)
h = @(theta) snc_env_poisson(lambda, theta);
end

function h = makeMap(D0, D1)
h = @(theta) snc_env_map(D0, D1, theta);
end

function h = makeSum(parts)
h = @(theta) sumEnv(theta, parts);
end

function h = makeLeftover(mu, cross)
h = @(theta) leftoverEnv(theta, mu, cross);
end

function h = makeOutput(arv, srv)
h = @(theta) outputEnv(theta, arv, srv);
end

function [sigma, rho] = sumEnv(theta, parts)
% Superposition of independent flows: the exponential forms multiply, so the
% envelope terms add.
sigma = 0; rho = 0;
for k = 1:numel(parts)
    [s, r] = parts{k}(theta);
    sigma = sigma + s;
    rho = rho + r;
end
end

function [sigma, rho] = leftoverEnv(theta, mu, cross)
[sS, rS] = snc_srv_exp(mu, theta);
if isempty(cross)
    sigma = sS; rho = rS;
    return
end
[sX, rX] = sumEnv(theta, cross);
[sigma, rho] = snc_leftover(sS, rS, sX, rX);
end

function [sigma, rho] = outputEnv(theta, arv, srv)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
[sigma, rho] = snc_output(sA, rA, sS, rS, theta);
end

function order = topoOrder(adj)
% Kahn's algorithm. Returns [] when the graph has a cycle, which is the
% feed-forward refusal.
n = size(adj,1);
indeg = sum(adj,1);
order = zeros(1,n);
placed = 0;
done = false(1,n);
while true
    cand = find(~done & indeg == 0, 1);
    if isempty(cand)
        break
    end
    placed = placed + 1;
    order(placed) = cand;
    done(cand) = true;
    indeg = indeg - double(adj(cand,:));
    indeg(cand) = 1; % keep it out of the candidate set
end
if placed < n
    order = [];
else
    order = order(1:placed);
end
end
