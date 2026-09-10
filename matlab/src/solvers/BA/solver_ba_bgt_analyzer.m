function [Q,U,R,T,C,X,lG,runtime] = solver_ba_bgt_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_BA_BGT_ANALYZER(SN, OPTIONS)
%
% Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
% multitype open Markovian network, valid for EVERY work-conserving Markovian
% policy. The polyhedron and the bound are NPFQN_BND_BGT; this analyzer maps
% the LINE model onto them and reads the bound back per station and class.
%
% CLASS SPACE. The reference's network is a MULTITYPE one: each type follows a
% FIXED sequence of stages, and stage k of type i is its own buffer. LINE's
% (station, job class) pair is that buffer, so the analyzer walks the routing
% matrix from the Source and turns each open class into one type whose stages
% are the pairs it visits. Two gates follow from the model and are enforced by
% name rather than approximated:
%
%   - ROUTING MUST BE DETERMINISTIC: a pair sends everything to one successor,
%     or everything to the Sink. A probabilistic split is a different network.
%   - ROUTES MUST NOT MERGE: a pair belongs to exactly one type. Where two
%     types share a buffer the reference's class index (i,k) is not defined,
%     and its arrival term L^j(i,1) lambda_i would be ambiguous.
%
% A re-entrant line is expressible: give the revisits distinct LINE classes
% (class switching), so each visit is its own pair.
%
% BOUND CONVENTION. Q(i,r) is the Theorem 4 bound; R follows by Little's law
% from it and the EXACT throughput T (an open network's per-class rates are
% fixed by the traffic equations, not by the policy), as does C. U is exact for
% the same reason.
%
% THE BOUND IS LOOSE, and knowingly so: the exception parameter of the smoothed
% Lyapunov function carries (Lmax+gamma)^3/gamma^2 and dominates as soon as
% there is more than one station. On M/M/1 it is 18x to 56x the exact mean
% queue length (tighter as the load rises); on a two-station tandem it is three
% orders of magnitude above. What is sharp is the STABILITY CERTIFICATE -- a
% feasible gamma > 0 proves every work-conserving policy stable, and the LP
% correctly refuses the Lu-Kumar network at per-station loads of 0.7, where
% global stability genuinely fails -- and the geometric tail RATE.
%
% Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance of
% multiclass Markovian queueing networks via piecewise linear Lyapunov
% functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1.
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
% run raises the sentence the report gave (the deterministic, non-merging route
% rule included); the checks below and in the walk are its guards.
reason = ba_open_refusal(sn, 'bgt.upper');
if ~isempty(reason)
    line_error(mfilename, '%s', reason);
end
if any(isfinite(sn.njobs))
    line_error(mfilename, ...
        'Method ''bgt.upper'' supports fully open networks only (no closed classes).');
end
isSource = false(M,1);
for i = 1:M
    isSource(i) = (sn.nodetype(sn.stationToNode(i)) == NodeType.Source);
end
if ~any(isSource)
    line_error(mfilename, 'Method ''bgt.upper'' requires an open network with a Source station.');
end
qstat = find(~isSource);
if any(sn.sched(qstat) == SchedStrategy.INF)
    line_error(mfilename, ...
        'Method ''bgt.upper'' does not support delay (infinite-server) stations: the reference''s network has one server per station.');
end
if any(sn.nservers(qstat) > 1)
    line_error(mfilename, 'Method ''bgt.upper'' does not support multi-server stations.');
end

rtst = sn_rt_stations(sn);
srcList = find(isSource);

% Pair (station, class) -> flat index in station space.
nq = numel(qstat);
pairStation = zeros(nq*K,1); pairClass = zeros(nq*K,1); pairFlat = zeros(nq*K,1);
np = 0;
for a = 1:nq
    i = qstat(a);
    for r = 1:K
        np = np + 1;
        pairStation(np) = i; pairClass(np) = r; pairFlat(np) = (i-1)*K + r;
    end
end
flatToPair = zeros(M*K,1);
flatToPair(pairFlat) = 1:np;

% ----- walk one deterministic route per source class -----
lambda = []; muCell = {}; sigmaCell = {}; routeCell = {}; srcClass = [];
used = false(np,1);
for si = 1:numel(srcList)
    s = srcList(si);
    for r0 = 1:K
        arr = sn.rates(s,r0);
        if ~isfinite(arr) || arr <= 0
            continue
        end
        row = full(rtst((s-1)*K + r0, :));
        cur = bgt_single_successor(row, flatToPair, pairFlat, ...
            sprintf('the Source for class %d', r0));
        route = [];
        while cur > 0
            if used(cur)
                line_error(mfilename, ...
                    ['Method ''bgt.upper'' needs routes that do not merge: station %d class %d ' ...
                     'is visited by more than one type. Give the visits distinct job classes.'], ...
                    pairStation(cur), pairClass(cur));
            end
            used(cur) = true;
            route(end+1) = cur; %#ok<AGROW>
            row = full(rtst(pairFlat(cur), :));
            cur = bgt_single_successor(row, flatToPair, pairFlat, ...
                sprintf('station %d class %d', pairStation(route(end)), pairClass(route(end))));
        end
        if isempty(route)
            line_error(mfilename, 'Class %d leaves the Source and reaches no station.', r0);
        end
        lambda(end+1,1) = arr; %#ok<AGROW>
        srcClass(end+1,1) = r0; %#ok<AGROW>
        routeCell{end+1} = route; %#ok<AGROW>
        muv = zeros(1,numel(route)); stv = zeros(1,numel(route));
        for k = 1:numel(route)
            p = route(k);
            muv(k) = sn.rates(pairStation(p), pairClass(p));
            stv(k) = pairStation(p);
            if ~isfinite(muv(k)) || muv(k) <= 0
                line_error(mfilename, ...
                    'Station %d has no service rate for class %d but carries its traffic.', ...
                    pairStation(p), pairClass(p));
            end
            if sn.procid(pairStation(p), pairClass(p)) ~= ProcessType.EXP
                line_error(mfilename, ...
                    'Method ''bgt.upper'' requires exponential service: station %d class %d is %s.', ...
                    pairStation(p), pairClass(p), ...
                    ProcessType.toText(sn.procid(pairStation(p), pairClass(p))));
            end
        end
        muCell{end+1} = muv; %#ok<AGROW>
        sigmaCell{end+1} = stv; %#ok<AGROW>
    end
end
if isempty(lambda)
    line_error(mfilename, 'The model carries no open traffic.');
end

% Dense station index space for the LP.
visited = unique(cell2mat(cellfun(@(v) v(:)', sigmaCell, 'UniformOutput', false)));
remap = zeros(M,1);
remap(visited) = 1:numel(visited);
for i = 1:numel(sigmaCell)
    sigmaCell{i} = remap(sigmaCell{i})';
    sigmaCell{i} = sigmaCell{i}(:)';
end

[Qub, info] = npfqn_bnd_bgt(lambda, muCell, sigmaCell, numel(visited));

% ----- read the bound back per station and class -----
for i = 1:numel(routeCell)
    route = routeCell{i};
    for k = 1:numel(route)
        p = route(k);
        ist = pairStation(p); r = pairClass(p);
        Q(ist,r) = Q(ist,r) + Qub{i}(k);
        T(ist,r) = T(ist,r) + lambda(i);
        U(ist,r) = U(ist,r) + lambda(i) / sn.rates(ist,r);
    end
end
for i = 1:M
    for r = 1:K
        if T(i,r) > 0
            R(i,r) = Q(i,r) / T(i,r);
        end
    end
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
for r = 1:K
    if X(1,r) > 0
        C(1,r) = sum(Q(:,r)) / X(1,r);
    end
end

line_debug('BGT bound: gamma=%g Lmax=%g B=%g U=%g tailRatio=%g', ...
    info.gamma, info.Lmax, info.B, info.U, info.tailRatio);
runtime = toc(T0);
end

% The single successor of a routing row, as a pair index, or 0 when everything
% leaves the network. A probabilistic split is refused by name: the reference's
% network has deterministic routing and a split is a different model, not an
% approximation of this one.
function p = bgt_single_successor(row, flatToPair, pairFlat, who)
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
    p = 0;   % everything leaves for the Sink
    return
end
if abs(mass - 1) > tol || abs(best - 1) > tol
    line_error('solver_ba_bgt_analyzer', ...
        ['Method ''bgt.upper'' needs deterministic routing: %s splits its departures ' ...
         '(the largest branch carries %.6g of them). The reference''s network routes ' ...
         'each type along a fixed sequence of stages.'], who, best);
end
end
