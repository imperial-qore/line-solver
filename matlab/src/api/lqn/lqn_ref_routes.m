%{
%{
 % @file lqn_ref_routes.m
 % @brief Synchronous call DAG carrying reference-task customers into a layer.
%}
%}

function [R, why] = lqn_ref_routes(lqn, callers, maxpaths, serverSet)
%{
%{
 % @brief Synchronous call DAG carrying reference-task customers into a layer.
 % @fn lqn_ref_routes(lqn, callers, maxpaths, serverSet)
 % @param lqn LayeredNetworkStruct.
 % @param callers Task indices that call the layer's server.
 % @param maxpaths (Optional) Refuse the layer above this many reference-task
 %        routes into it. Default 32. The count is derived, never enumerated.
 % @param serverSet (Optional) Server elements of the layer. A prefix node whose
 %        task is one of them is recursion and refuses the layer.
 % @return R Struct array, one element per reference task reaching CALLERS.
 % @return WHY Non-empty when the layer must fall back to another interlock
 %        method; R is then empty.
%}
%}
% [R, WHY] = LQN_REF_ROUTES(LQN, CALLERS, MAXPATHS, SERVERSET)
%
% Resolves, for the caller set of one layer, the synchronous call graph along
% which reference (REF) task customers descend to those callers, and the mean
% number of times each entry and each call is invoked per REF cycle.
%
% SolverLN gives every caller of a layer its own client chain, each carrying its
% own thread pool as a population, so a layer holds one customer per caller even
% when those callers are the same REF customers arriving by different routes.
% The caller set of this group, .members, is what has to become ONE chain.
%
% The result carries a single DAG and two readings of it:
%
%   .members    every caller of the layer on the DAG, in descent order. The
%               first is the chain head when the REF task is itself a caller.
%   .prefixPos  the entries from the REF task down to, but not past, the first
%               caller on each path. Below that first caller the layer already
%               holds an activity subgraph, and the descent through it is
%               spliced call by call, so the prefix has nothing left to say.
%
% Nodes are ENTRIES, because the work a hop charges depends on which entry was
% called. Visits are computed TOPOLOGICALLY, v(u) = sum over parents of
% v(p)*w(a)*callmean, and routes are COUNTED in the same pass rather than
% enumerated: the weights a route decomposition would produce are obtainable in
% one linear pass, and two sources of truth for one quantity is one too many.
%
% Only SYNC calls are followed. An ASYNC call is send-no-reply, so it terminates
% blocking and the customer below it is not the REF customer; forwarding is
% already flattened into pseudo-SYNC arcs by lqn_fwd_rendezvous before the
% layers are built, and any CallType.FWD left in the struct carries no blocking.
%
% R(g) carries:
%   .reftask       task index of the REF task at the root
%   .headIsCaller  true when the REF task is a caller of this layer
%   .members       callers of this layer on the DAG, descent order
%   .entries       every DAG entry, topologically ordered, root entries first
%   .etask         lqn.parent of each entry
%   .ismember      true where that entry's task is a caller of this layer
%   .vEntry        mean invocations of each entry per REF cycle
%   .actweight     cell, per entry, [aidx; executions per invocation]
%   .calls         [cidx, fromPos, toPos, aidx, vCall] rows, vCall being the
%                  mean invocations of that call per REF cycle
%   .prefixPos     positions in .entries forming the prefix, topological order
%   .prefixTerm    true where that prefix position is a first caller
%   .npaths        distinct REF-to-caller routes, counted
%   .poolmin       min of lqn.maxmult over the DAG tasks. A DIAGNOSTIC: the
%                  chain population is never capped by it.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = struct('reftask',{},'headIsCaller',{},'members',{},'entries',{},'etask',{}, ...
    'ismember',{},'vEntry',{},'actweight',{},'calls',{},'prefixPos',{}, ...
    'prefixTerm',{},'npaths',{},'poolmin',{});
why = '';

if nargin < 3 || isempty(maxpaths)
    maxpaths = 32;
end
if nargin < 4
    serverSet = [];
end
if ~isfield(lqn,'ncalls') || lqn.ncalls == 0 || isempty(callers)
    return
end

callers = unique(callers(:))';
isCaller = false(lqn.nidx,1);
isCaller(callers) = true;

succ = syncSuccessors(lqn);

reftasks = [];
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if full(lqn.isref(tidx))
        reftasks(end+1) = tidx; %#ok<AGROW>
    end
end
if isempty(reftasks)
    return
end

% Which REF tasks reach which callers. A caller reachable from two REF tasks is
% two INDEPENDENT customer pools, and merging them into one chain would invent a
% correlation the model does not contain. The whole LAYER falls back rather than
% that caller alone: a refused caller may simultaneously lie on another group's
% path, which would give the same threads both a fallback chain of their own and
% a spliced representation, i.e. count them twice.
nrefOf = zeros(lqn.nidx,1);
for g = 1:numel(reftasks)
    seen = reachableEntries(lqn, succ, reftasks(g));
    mem = unique(full(lqn.parent(find(seen)))); %#ok<FNDSB>
    mem = mem(isCaller(mem));
    if isCaller(reftasks(g))
        mem = union(mem, reftasks(g));
    end
    nrefOf(mem) = nrefOf(mem) + 1;
end
overloaded = find(nrefOf > 1)';
if ~isempty(overloaded)
    why = sprintf(['task ''%s'' is reachable from %d reference tasks, whose customer ' ...
        'pools are independent'], nameOf(lqn, overloaded(1)), full(nrefOf(overloaded(1))));
    return
end

for g = 1:numel(reftasks)
    [grp, gwhy] = buildGroup(lqn, succ, reftasks(g), isCaller, maxpaths, serverSet);
    if ~isempty(gwhy)
        why = gwhy;
        R = R([]);
        return
    end
    if ~isempty(grp)
        R(end+1) = grp; %#ok<AGROW>
    end
end
end

function succ = syncSuccessors(lqn)
% Per calling ENTRY, the rows [cidx, called entry, calling activity, callmean]
% of its synchronous calls. callpair(:,1) is the activity issuing the call and
% callpair(:,2) the entry called.
%
% The calling entry is NOT lqn.parent of the activity: an activity's parent is
% its TASK, not its entry (in 10-interlock, parent(A:e0_ph1) is R:t0, not E:e0).
% Keying the successor lists on the parent files every call under a task index
% while the traversal reads entry indices, and every sweep comes back empty.
% Invert actsof over the ENTRY range instead.
succ = cell(lqn.nidx,1);
entryOfAct = zeros(lqn.nidx,1);
for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    acts = lqn.actsof{eidx};
    if ~isempty(acts)
        entryOfAct(acts) = eidx;
    end
end
for cidx = 1:lqn.ncalls
    if full(lqn.calltype(cidx)) ~= CallType.SYNC
        continue
    end
    aidx = full(lqn.callpair(cidx,1));
    eidx_to = full(lqn.callpair(cidx,2));
    if aidx < 1 || eidx_to < 1
        continue
    end
    eidx_from = entryOfAct(aidx);
    if eidx_from < 1
        continue
    end
    succ{eidx_from}(end+1,:) = [cidx, eidx_to, aidx, callMeanOf(lqn, cidx)];
end
end

function m = callMeanOf(lqn, cidx)
% Mean number of invocations carried by call CIDX
m = NaN;
if isfield(lqn,'callproc_mean') && numel(lqn.callproc_mean) >= cidx
    m = full(lqn.callproc_mean(cidx));
end
if ~isfinite(m) && isfield(lqn,'callproc') && numel(lqn.callproc) >= cidx && ~isempty(lqn.callproc{cidx})
    m = lqn.callproc{cidx}.getMean;
end
if ~isfinite(m)
    m = 1;
end
end

function nm = nameOf(lqn, idx)
% Printable name of an LQN element
if isfield(lqn,'hashnames') && numel(lqn.hashnames) >= idx && ~isempty(lqn.hashnames{idx})
    nm = lqn.hashnames{idx};
else
    nm = sprintf('#%d', idx);
end
end

function seen = reachableEntries(lqn, succ, tidx)
% Every entry reachable from task TIDX over SYNC calls, without pruning: a
% caller of the layer may lie below another caller, and both belong to the
% merged chain.
seen = false(lqn.nidx,1);
stack = lqn.entriesof{tidx}(:)';
while ~isempty(stack)
    eidx = stack(end);
    stack(end) = [];
    if seen(eidx)
        continue
    end
    seen(eidx) = true;
    for k = 1:size(succ{eidx},1)
        stack(end+1) = succ{eidx}(k,2); %#ok<AGROW>
    end
end
end

function [g, why] = buildGroup(lqn, succ, r, isCaller, maxpaths, serverSet)
% One group: the full synchronous DAG below REF task R, its per-cycle visit
% counts, and the prefix that sits above the first callers of the layer
g = [];
why = '';

roots = lqn.entriesof{r}(:)';
if isempty(roots)
    return
end

% Depth-first sweep with an on-stack marker. A back edge is REFUSED rather than
% pruned the way traceInterlockPaths prunes for its flow table: a cycle makes
% v(u) a geometric series in a quantity the chain would have to express as a
% self-loop through the layer's own server.
WHITE = 0; GREY = 1; BLACK = 2;
color = zeros(lqn.nidx,1);
post = [];
for e0 = roots
    if color(e0) ~= WHITE
        continue
    end
    stack = [e0, 0];
    while ~isempty(stack)
        u = stack(end,1);
        ci = stack(end,2);
        if ci == 0
            color(u) = GREY;
        end
        kids = zeros(1,0);
        if ~isempty(succ{u})
            kids = succ{u}(:,2)';
        end
        if ci < numel(kids)
            stack(end,2) = ci + 1;
            v = kids(ci+1);
            if color(v) == GREY
                why = sprintf('the synchronous call graph below ''%s'' is cyclic', nameOf(lqn, v));
                return
            elseif color(v) == WHITE
                stack(end+1,:) = [v, 0]; %#ok<AGROW>
            end
        else
            color(u) = BLACK;
            post(end+1) = u; %#ok<AGROW>
            stack(end,:) = [];
        end
    end
end

entries = fliplr(post);              % topological order
if isempty(entries)
    return
end
pos = zeros(lqn.nidx,1);
pos(entries) = 1:numel(entries);
etask = full(lqn.parent(entries));
etask = etask(:)';
ismem = isCaller(etask)';

if ~any(ismem)
    return % this REF reaches none of the callers
end

% Per-entry activity weights, then the call list
actweight = cell(1, numel(entries));
calls = zeros(0,5);
for i = 1:numel(entries)
    u = entries(i);
    [aw, awwhy] = entryActWeights(lqn, u);
    if ~isempty(awwhy)
        why = awwhy;
        return
    end
    actweight{i} = aw;
    for k = 1:size(succ{u},1)
        v = succ{u}(k,2);
        if pos(v) == 0
            continue
        end
        aidx = succ{u}(k,3);
        w = 0;
        if ~isempty(aw)
            hit = find(aw(1,:) == aidx, 1);
            if ~isempty(hit)
                w = aw(2, hit);
            end
        end
        calls(end+1,:) = [succ{u}(k,1), i, pos(v), aidx, w * succ{u}(k,4)]; %#ok<AGROW>
    end
end

% Topological visits. The REF task selects among its own entries the way
% buildLayersRecursive already splits a caller class, with equal probability
% refreshed to the throughput ratio, so a root entry is visited 1/nentries
% times per REF cycle.
vEntry = zeros(1, numel(entries));
rootpos = pos(roots);
rootpos = rootpos(rootpos > 0)';
vEntry(rootpos) = 1 / numel(roots);
for i = 1:numel(entries)
    for k = find(calls(:,2)' == i)
        vEntry(calls(k,3)) = vEntry(calls(k,3)) + vEntry(i) * calls(k,5);
    end
end
vCall = zeros(size(calls,1),1);
for k = 1:size(calls,1)
    vCall(k) = vEntry(calls(k,2)) * calls(k,5);
end
calls(:,5) = vCall;

% The prefix: entries reachable from a root without passing THROUGH a caller.
% A caller entry reached this way is the prefix's terminal, where the layer's
% own activity subgraph takes over.
inPrefix = false(1, numel(entries));
prefTerm = false(1, numel(entries));
npath = zeros(1, numel(entries));
inPrefix(rootpos) = true;
npath(rootpos) = 1;
for i = 1:numel(entries)
    if ~inPrefix(i)
        continue
    end
    if ismem(i)
        prefTerm(i) = true;
        continue % do not descend past a caller
    end
    % A HOP whose task is a server of this layer would charge at the client
    % Delay work that the layer's own station exists to serve, and the customer
    % would be in two places at once. A node BELOW the first caller is not a hop
    % (the layer's activity subgraph serves it) and a terminal is not a hop
    % either, so only this position is recursion.
    if ~isempty(serverSet) && any(serverSet == etask(i))
        why = sprintf(['task ''%s'' is both an intermediate on the reference path and ' ...
            'a server of this layer'], nameOf(lqn, etask(i)));
        return
    end
    for k = find(calls(:,2)' == i)
        j = calls(k,3);
        inPrefix(j) = true;
        npath(j) = npath(j) + npath(i);
    end
end

np = sum(npath(prefTerm));
if np > maxpaths
    why = sprintf(['the reference path into this layer carries %g distinct routes, ' ...
        'above config.interlock_maxpaths = %g'], np, maxpaths);
    return
end

g = struct( ...
    'reftask', r, ...
    'headIsCaller', isCaller(r), ...
    'members', unique(etask(ismem), 'stable'), ...
    'entries', entries, ...
    'etask', etask, ...
    'ismember', ismem, ...
    'vEntry', vEntry, ...
    'actweight', {actweight}, ...
    'calls', calls, ...
    'prefixPos', find(inPrefix), ...
    'prefixTerm', prefTerm(inPrefix), ...
    'npaths', np, ...
    'poolmin', poolMin(lqn, etask));
end

function p = poolMin(lqn, tasks)
% Smallest thread pool on the path. DIAGNOSTIC ONLY: capping the chain at this
% would delete customers from the REF task's own think stage, not merely from
% the in-flight portion, and would understate throughput in the direction
% opposite to the error interlocking exists to remove.
p = Inf;
if ~isfield(lqn,'maxmult')
    return
end
for t = unique(tasks(tasks > 0))
    if numel(lqn.maxmult) >= t
        p = min(p, full(lqn.maxmult(t)));
    end
end
end

function [aw, why] = entryActWeights(lqn, eidx)
% AW is [aidx; weight], the mean executions of each activity of entry EIDX per
% invocation of that entry, propagated over lqn.graph so an OR-branch splits its
% successors by the declared probabilities instead of counting both in full.
aw = zeros(2,0);
why = '';
acts = lqn.actsof{eidx};
if isempty(acts)
    return
end
acts = acts(:)';
nodeset = [eidx, acts];
pos = zeros(lqn.nidx,1);
pos(nodeset) = 1:numel(nodeset);

A = zeros(numel(nodeset));
for i = 1:numel(nodeset)
    u = nodeset(i);
    for v = find(lqn.graph(u,:))
        if pos(v) > 0
            A(i, pos(v)) = full(lqn.graph(u,v));
        end
    end
end

% Topological propagation with an explicit cycle test. An activity-graph loop
% makes the executions a geometric series the chain cannot express, so the layer
% falls back rather than silently truncating the count.
indeg = sum(A > 0, 1);
indeg(1) = 0; % the entry is the source
w = zeros(1, numel(nodeset));
w(1) = 1;
remaining = indeg;
queue = find(remaining == 0);
done = false(1, numel(nodeset));
ndone = 0;
while ~isempty(queue)
    i = queue(1);
    queue(1) = [];
    if done(i)
        continue
    end
    done(i) = true;
    ndone = ndone + 1;
    for j = find(A(i,:) > 0)
        w(j) = w(j) + w(i) * A(i,j);
        remaining(j) = remaining(j) - 1;
        if remaining(j) <= 0 && ~done(j)
            queue(end+1) = j; %#ok<AGROW>
        end
    end
end
if ndone < numel(nodeset)
    why = sprintf('the activity graph of entry ''%s'' contains a loop', nameOf(lqn, eidx));
    return
end

aw = [acts; w(2:end)];
end
