function reason = ba_open_refusal(sn, method)
% REASON = BA_OPEN_REFUSAL(SN, METHOD)
%
% The structural premises of the three OPEN-network bound families, 'bpt',
% 'bgt' and 'snc', in one place: the reason the resolved METHOD cannot bound
% SN, or '' when it can. Empty for every other method.
%
% ONE PREDICATE, TWO CALLERS. Each of solver_ba_bpt_analyzer,
% solver_ba_bgt_analyzer and solver_ba_snc_analyzer asks it before touching the
% model and raises what it returns; BA_METHOD_REFUSAL asks it so that
% supportsModelMethod, findSolver and listValidMethods report the same
% sentence. Until it existed the routing rules lived in the analyzers alone, so
% model.help called 'bgt.upper' runnable on an open model with a probabilistic
% split and 'snc.upper' runnable on a station graph with a cycle, and both
% raised the moment they were run.
%
% WHAT IS CHECKED, family by family. All three: a fully open model with a
% Source, one server per queueing station and no delay (each bound is derived
% for one exponential server per station), and exponential service at every
% (station, class) pair that CARRIES TRAFFIC -- the traffic equations are
% solved here, as the analyzers solve them, so a service time given to a pair
% its class never visits is not held against the model. 'bgt' additionally
% needs DETERMINISTIC, NON-MERGING routes: a pair sends everything to one
% successor or everything to the Sink, and belongs to exactly one type. 'snc'
% needs deterministic routing downstream of the Source (a split AT a Poisson
% Source is exact and allowed, a split of any other Markovian source is not),
% one service rate per station across the classes it serves, and a
% FEED-FORWARD station graph. The arrival law is deliberately NOT checked: 'snc'
% consumes any Markovian (D0,D1) source, and for 'bpt'/'bgt' the exponential
% arrival premise is a registry-expressible delta in SolverBA.getMethodFeatureSet.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
dot = strfind(method,'.');
if isempty(dot)
    fam = method;
else
    fam = method(1:dot(1)-1);
end
if ~any(strcmp(fam, {'bpt','bgt','snc'}))
    return
end

M = sn.nstations;
K = sn.nclasses;
tol = 1e-10;

if any(isfinite(sn.njobs))
    reason = sprintf('Method ''%s'' supports fully open networks only (no closed classes).', method);
    return
end
isSource = false(M,1);
for i = 1:M
    isSource(i) = (sn.nodetype(sn.stationToNode(i)) == NodeType.Source);
end
if ~any(isSource)
    reason = sprintf('Method ''%s'' requires an open network with a Source station.', method);
    return
end
qstat = find(~isSource);
if any(sn.sched(qstat) == SchedStrategy.INF)
    reason = sprintf(['Method ''%s'' does not support delay (infinite-server) stations: the ' ...
        'bound is derived for one server per station.'], method);
    return
end
if any(sn.nservers(qstat) > 1)
    reason = sprintf('Method ''%s'' does not support multi-server stations.', method);
    return
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
srcList = find(isSource);
inject = zeros(np, numel(srcList)*K);
lambda0 = zeros(np,1);
srcOfCol = zeros(1, numel(srcList)*K);
clsOfCol = zeros(1, numel(srcList)*K);
col = 0;
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
lamAll = (eye(np) - P') \ lambda0;
keep = find(lamAll > 1e-12 * max(1, max(lamAll)));
if isempty(keep)
    reason = 'The model carries no open traffic.';
    return
end

% ----- exponential service at the pairs that carry traffic -----
for kk = 1:numel(keep)
    p = keep(kk);
    i = pairStation(p); r = pairClass(p);
    mu = sn.rates(i,r);
    if ~isfinite(mu) || mu <= 0
        reason = sprintf('Station %d has no service rate for class %d but carries its traffic.', i, r);
        return
    end
    if sn.procid(i,r) ~= ProcessType.EXP
        reason = sprintf('Method ''%s'' requires exponential service: station %d class %d is %s.', ...
            method, i, r, ProcessType.toText(sn.procid(i,r)));
        return
    end
end

switch fam
    case 'bgt'
        % One deterministic route per (source, class), and no pair on two routes.
        used = false(np,1);
        for si = 1:numel(srcList)
            s = srcList(si);
            for r0 = 1:K
                arr = sn.rates(s,r0);
                if ~isfinite(arr) || arr <= 0
                    continue
                end
                row = full(rtst((s-1)*K + r0, :));
                [cur, splitWhy] = bgt_successor(row, pairFlat, sprintf('the Source for class %d', r0));
                if ~isempty(splitWhy)
                    reason = sprintf('Method ''%s'' needs deterministic routing: %s', method, splitWhy);
                    return
                end
                if cur == 0
                    reason = sprintf('Class %d leaves the Source and reaches no station.', r0);
                    return
                end
                while cur > 0
                    if used(cur)
                        reason = sprintf(['Method ''%s'' needs routes that do not merge: station %d ' ...
                            'class %d is visited by more than one type. Give the visits distinct job ' ...
                            'classes.'], method, pairStation(cur), pairClass(cur));
                        return
                    end
                    used(cur) = true;
                    row = full(rtst(pairFlat(cur), :));
                    [nxt, splitWhy] = bgt_successor(row, pairFlat, ...
                        sprintf('station %d class %d', pairStation(cur), pairClass(cur)));
                    if ~isempty(splitWhy)
                        reason = sprintf('Method ''%s'' needs deterministic routing: %s', method, splitWhy);
                        return
                    end
                    cur = nxt;
                end
            end
        end
    case 'snc'
        Pk = P(keep,keep);
        injk = inject(keep,:);
        stK = pairStation(keep);
        clK = pairClass(keep);
        nk = numel(keep);
        for p = 1:nk
            succ = find(Pk(p,:) > tol);
            if numel(succ) > 1
                reason = sprintf(['Method ''%s'' requires deterministic routing downstream of the ' ...
                    'Source: station %d class %d splits its flow over %d destinations.'], ...
                    method, stK(p), clK(p), numel(succ));
                return
            end
            if numel(succ) == 1 && abs(Pk(p,succ) - 1) > 1e-8
                reason = sprintf(['Method ''%s'' requires deterministic routing downstream of the ' ...
                    'Source: station %d class %d routes onward with probability %g.'], ...
                    method, stK(p), clK(p), Pk(p,succ));
                return
            end
        end
        for c = 1:size(injk,2)
            s = srcOfCol(c); r0 = clsOfCol(c);
            if s == 0 || sum(injk(:,c)) <= 0
                continue
            end
            dest = find(injk(:,c) > tol);
            if numel(dest) > 1 && sn.procid(s,r0) ~= ProcessType.EXP
                reason = sprintf(['Method ''%s'' can split only a Poisson Source: source %d class %d ' ...
                    'is %s and feeds %d stations.'], method, s, r0, ...
                    ProcessType.toText(sn.procid(s,r0)), numel(dest));
                return
            end
        end
        stationsUsed = unique(stK);
        for a = 1:numel(stationsUsed)
            i = stationsUsed(a);
            sel = find(stK == i);
            rates = zeros(numel(sel),1);
            for b = 1:numel(sel)
                rates(b) = sn.rates(i, clK(sel(b)));
            end
            if max(rates) - min(rates) > 1e-8 * max(1, max(rates))
                reason = sprintf(['Method ''%s'' requires the classes sharing a station to have equal ' ...
                    'service rates: station %d carries rates in [%g, %g].'], method, i, min(rates), max(rates));
                return
            end
        end
        ns = numel(stationsUsed);
        stIdx = zeros(M,1);
        stIdx(stationsUsed) = 1:ns;
        adj = false(ns,ns);
        for p = 1:nk
            succ = find(Pk(p,:) > tol);
            if isempty(succ)
                continue
            end
            adj(stIdx(stK(p)), stIdx(stK(succ))) = true;
        end
        if ~isAcyclic(adj)
            reason = sprintf(['Method ''%s'' requires a feed-forward network: the station graph ' ...
                'has a cycle, so a station''s cross traffic is not determined upstream of it.'], method);
            return
        end
end
end

function [p, why] = bgt_successor(row, pairFlat, who)
% The single successor of a routing row as a pair index, 0 when everything
% leaves for the Sink, and WHY when the row splits. Same tolerance and same
% test as bgt_single_successor in solver_ba_bgt_analyzer.
why = '';
tol = 1e-9;
mass = 0; p = 0; best = 0;
for q = 1:numel(pairFlat)
    v = row(pairFlat(q));
    if v > tol
        mass = mass + v;
        if v > best
            best = v; p = q;
        end
    end
end
if mass <= tol
    p = 0;
    return
end
if abs(mass - 1) > tol || abs(best - 1) > tol
    why = sprintf(['%s splits its departures (the largest branch carries %.6g of them). The ' ...
        'reference''s network routes each type along a fixed sequence of stages.'], who, best);
end
end

function bool = isAcyclic(adj)
% Kahn's algorithm, as topoOrder in solver_ba_snc_analyzer: true when every
% node can be placed, false when a cycle stops the placement.
n = size(adj,1);
indeg = sum(adj,1);
done = false(1,n);
placed = 0;
while true
    cand = find(~done & indeg == 0, 1);
    if isempty(cand)
        break
    end
    placed = placed + 1;
    done(cand) = true;
    indeg = indeg - double(adj(cand,:));
    indeg(cand) = 1;
end
bool = (placed == n);
end
