function model = LQN2QN(lqn)
% LQN2QN Convert a LayeredNetwork (LQN) to a Network (QN) using REPLY signals
%
% model = LQN2QN(lqn) flattens a LayeredNetwork into a single queueing
% network in which synchronous call blocking is represented by REPLY
% signals.
%
% Construction
% - One station per host processor (scheduling and multiplicity taken from
%   the processor). Tasks sharing a processor share the station, as in the
%   LQN semantics where the processor is the contended resource.
% - One Delay per reference task, holding its think time.
% - One closed class per step of the expanded activity graph. A step is an
%   activity, or one call stage of an activity that issues synchronous
%   calls. Steps of the reference task chain carry population 0 except the
%   think class, which carries the reference task multiplicity.
% - A synchronous call site blocks its caller: the step class has a REPLY
%   signal bound to it (sn.syncreply), the token proceeds to the callee,
%   and the callee's replying activity class-switches to that signal, which
%   returns to the caller station and unblocks it.
% - Call multiplicity mean m is unrolled into floor(m) mandatory call
%   stages plus, if m is not integer, one further stage entered with
%   probability m-floor(m). This preserves both the mean number of calls
%   and the blocking of each individual call.
% - OR-branch and loop probabilities are read from the activity graph
%   weights lsn.graph(a,b). Multiple entries per task and branching call
%   trees are expanded, each call site receiving its own copy of the callee
%   subgraph so that per-call-path response times remain distinguishable.
%
% - AND precedences become Fork and Join nodes: a POST_AND successor set is
%   entered through a Fork and one Router per branch, since a Fork cannot
%   switch class per output link, and the PRE_AND branch tails switch back
%   to the class that entered the Fork, which is the class the Join matches
%   its siblings on. A branch tail that issues a synchronous call is given a
%   merge step, so that it reaches the Join in an ordinary class rather than
%   as a REPLY signal, which carries no forked-task identity.
%
% - A CacheTask becomes a Cache node. The activity bound to an ItemEntry is
%   the read step and sits on that node, carrying the item popularity and
%   cardinality of the entry; its two CacheAccess successors become the hit
%   and the miss class, and the class switch is performed by the Cache node
%   itself, so the routes leaving it are written in the successor class.
%
% - An asynchronous call is lowered to a non-blocking visit: the caller does
%   not hold its server for the duration of the call, but it is serialised
%   behind it, since a closed network has no means of creating the second
%   token that a truly concurrent send would require.
%
% - Entry forwarding splits the reply exits of the forwarding entry: with
%   the forwarding probability the request is handed to the target entry,
%   which replies to the original caller, so the forwarder is released while
%   the caller stays blocked.
%
% - The multiplicity of a non-reference task is its thread pool: at most
%   that many requests may be inside the task at once, where inside spans
%   the task's own steps and those of its nested callees, since a thread is
%   held for the whole of a synchronous call. It is enforced by a finite
%   capacity region with one linear admission constraint per task, and the
%   calls of such a task do not hold the caller's server, since the
%   processor is released while a thread waits for a reply.
%
% - An entry with an open arrival process receives requests from a Source:
%   they traverse the entry's subgraph in open classes and leave through a
%   Sink where the entry would reply. Calls on an open chain do not hold
%   the caller's server, since a REPLY signal is a closed class that cannot
%   be woven into an open chain; the thread pools of the traversed tasks
%   are still enforced by the finite capacity region.
%
% - Phase-2 activities, the successors of a replying activity, run after the
%   reply: the replying step's exit routes back to the caller, and each of
%   its service completions spawns the continuation at the host station
%   (sn.classspawn). The spawned token walks the phase-2 subgraph holding
%   only the task's own thread and is destroyed at the chain end, through a
%   NEGATIVE signal that always misses on a closed chain, or through the
%   Sink on an open one. A boundary that ends on a call site is
%   normalised through a merge step at the host station; a phase 2 that
%   opens with an AND-fork spawns into an immediate head that feeds the
%   Fork; at an AND-join branch tail the spawned token inherits the fork
%   identity of the trigger and stands in for it at the Join; at a cache
%   read the reply is emitted by an immediate trigger step per hit/miss
%   outcome, whose completion spawns the matching branch continuation.
%
% Not yet represented: delayed-hit retrieval on the cache miss path, and the
% thread pool of a task with an internal AND-fork. Each is reported through
% line_warning.
%
% Example:
%   lqn = LayeredNetwork('MyLQN');
%   % ... define LQN model ...
%   model = LQN2QN(lqn);
%   SolverLDES(model).getAvgTable()
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lsn = lqn.getStruct();
model = Network([lqn.getName(), '-QN']);

refTaskIndices = find(lsn.isref);

%% Entries with an open arrival process
% Source -> open classes -> Sink; open-chain calls never hold the caller's server -- see _kb/04-networkstruct.md
openEntries = [];
if isfield(lsn, 'arrival') && ~isempty(lsn.arrival)
    for e_ = lsn.eshift+1:lsn.eshift+lsn.nentries
        if e_ <= length(lsn.arrival) && isa(lsn.arrival{e_}, 'Distribution') && ...
                isfinite(lsn.arrival{e_}.getMean()) && lsn.arrival{e_}.getMean() > GlobalConstants.FineTol
            openEntries(end+1) = e_; %#ok<AGROW>
        end
    end
end

if isempty(refTaskIndices) && isempty(openEntries)
    line_error(mfilename, 'LQN must have at least one reference task or open arrival.');
end

MAXCALLSTAGES = 20;   % guard against unrolling a huge call multiplicity

%% Unsupported features, reported once each
warnUnsupported(lsn);

%% Tasks whose multiplicity is a thread pool
% One thread per request, capped by a finite capacity region instead of server-hold blocking -- see _kb/04-networkstruct.md
fcrTask = false(lsn.nhosts + lsn.ntasks, 1);
for t_ = 1:lsn.ntasks
    tidx_ = lsn.tshift + t_;
    if lsn.isref(tidx_) || (tidx_ <= length(lsn.iscache) && lsn.iscache(tidx_))
        continue;
    end
    if ~isfinite(lsn.mult(tidx_)) || lsn.sched(tidx_) == SchedStrategy.INF
        continue;
    end
    if taskHasAndFork(tidx_)
        line_warning(mfilename, sprintf(['Multiplicity of task %s is not enforced: ' ...
            'an AND-fork inside a task cannot be capped by a finite capacity ' ...
            'region, whose job count would double-count the forked siblings.'], ...
            lsn.names{tidx_}));
        continue;
    end
    fcrTask(tidx_) = true;
end

% A task called transitively from inside an AND-fork branch is excluded too -- see _kb/04-networkstruct.md
branchActs = [];
if isfield(lsn, 'actposttype') && ~isempty(lsn.actposttype)
    for a_ = lsn.ashift+1:lsn.ashift+lsn.nacts
        if a_ > length(lsn.actposttype) || full(lsn.actposttype(a_)) ~= ActivityPrecedenceType.POST_AND
            continue;
        end
        frontier_ = a_;
        while ~isempty(frontier_)
            cur_ = frontier_(1); frontier_(1) = [];
            if any(branchActs == cur_)
                continue;
            end
            branchActs(end+1) = cur_; %#ok<AGROW>
            if isfield(lsn, 'actpretype') && ~isempty(lsn.actpretype) && ...
                    cur_ <= length(lsn.actpretype) && ...
                    full(lsn.actpretype(cur_)) == ActivityPrecedenceType.PRE_AND
                continue;   % branch tail: do not traverse past the join
            end
            succ_ = find(lsn.graph(cur_, :));
            succ_ = succ_(succ_ > lsn.ashift & lsn.parent(succ_).' == lsn.parent(cur_));
            frontier_ = [frontier_, succ_]; %#ok<AGROW>
        end
    end
end
if ~isempty(branchActs) && any(fcrTask)
    front_ = [];
    for a_ = branchActs
        if a_ <= length(lsn.callsof) && ~isempty(lsn.callsof{a_})
            for c_ = lsn.callsof{a_}
                front_(end+1) = lsn.parent(lsn.callpair(c_, 2)); %#ok<AGROW>
            end
        end
    end
    shadow_ = false(size(fcrTask));
    while ~isempty(front_)
        t_ = front_(1); front_(1) = [];
        if shadow_(t_)
            continue;
        end
        shadow_(t_) = true;
        for a_ = lsn.ashift+1:lsn.ashift+lsn.nacts
            if lsn.parent(a_) == t_ && a_ <= length(lsn.callsof) && ~isempty(lsn.callsof{a_})
                for c_ = lsn.callsof{a_}
                    front_(end+1) = lsn.parent(lsn.callpair(c_, 2)); %#ok<AGROW>
                end
            end
        end
    end
    for t_ = find(shadow_(:).' & fcrTask(:).')
        fcrTask(t_) = false;
        line_warning(mfilename, sprintf(['Multiplicity of task %s is not enforced: ' ...
            'it is called from inside an AND-fork branch, whose flows the ' ...
            'fork-join transformation retags outside the admission constraint.'], ...
            lsn.names{t_}));
    end
end

%% Stations: one per host processor
hostStation = cell(lsn.nhosts, 1);
hostIsDelay = false(lsn.nhosts, 1);
for h = 1:lsn.nhosts
    nservers = lsn.mult(h);
    sched = lsn.sched(h);
    if isinf(nservers) || sched == SchedStrategy.INF
        hostStation{h} = Delay(model, lsn.names{h});
        hostIsDelay(h) = true;
    else
        q = Queue(model, lsn.names{h}, sched);
        q.setNumberOfServers(nservers);
        hostStation{h} = q;
    end
end

%% Think delays, one per reference task
thinkNode = containers.Map('KeyType', 'double', 'ValueType', 'any');
for rt = 1:length(refTaskIndices)
    refTidx = refTaskIndices(rt);
    thinkNode(refTidx) = Delay(model, [lsn.names{refTidx}, '_Think']);
end

%% Pass 1: expand the activity graph into a step graph
% Steps carry no LINE objects yet, so all classes can be created before the routing matrix is initialised
stepAidx = [];        % activity index of the step
stepHost = [];        % host processor index of the step station
stepSvc = {};         % service distribution at the step station, [] if none
stepName = {};        % class name
stepBlocks = [];      % true if the step blocks on a synchronous call
stepIsThink = [];     % true for a think step (reference task)
stepRefTask = [];     % reference task the step belongs to
stepNode = {};        % Fork/Join/Router node the step sits on, [] for stations
stepClassOwner = [];  % step whose class this step travels in (itself, normally)
stepTasks = {};       % thread-pool tasks holding a thread while at this step

% flow/reply/spawnPairs/ph2Exits step-graph arrays -- see _kb/04-networkstruct.md
flow = zeros(0, 5);
reply = zeros(0, 4);
spawnPairs = zeros(0, 2);
ph2Exits = zeros(0, 4);

entryStack = [];   % guards against recursive call cycles
threadStack = [];  % tasks holding a thread during expansion; tracks entryStack except across forwarding (releases the forwarder's thread)

% Cache nodes, one per CacheTask, and the read/hit/miss wiring to apply once
% the classes exist.
cacheNodeOf = containers.Map('KeyType', 'double', 'ValueType', 'any');
cacheWiring = struct('node', {}, 'readStep', {}, 'hitStep', {}, 'missStep', {}, ...
    'itemproc', {}, 'nitems', {});

for rt = 1:length(refTaskIndices)
    refTidx = refTaskIndices(rt);
    thinkStep = addStep([], [], [], [lsn.names{refTidx}, '_Think'], false, true, refTidx);

    entries = lsn.entriesof{refTidx};
    for eidx = entries
        [firstStep, replyExits, terminals] = expandEntry(eidx, refTidx);
        if isempty(firstStep)
            continue;
        end
        addRoute([thinkStep, 0], firstStep, 1.0);
        % A reference task has no caller: its replies and its dead ends
        % both close the cycle at the think delay.
        exits = [replyExits; terminals];
        for s = 1:size(exits, 1)
            addRoute(exits(s, 1:2), thinkStep, exits(s, 3));
        end
    end
end

%% Pass 1b: open arrival chains
% openWiring rows: {eidx, firstStep, exits}, wired to Source/Sink in pass 4; reference task index 0 marks an open-chain step
openWiring = struct('eidx', {}, 'firstStep', {}, 'exits', {});
srcNode = []; snkNode = [];
if ~isempty(openEntries)
    srcNode = Source(model, 'Source');
    snkNode = Sink(model, 'Sink');
end
for eidx = openEntries
    [firstStep, replyExits, terminals] = expandEntry(eidx, 0);
    if isempty(firstStep)
        line_warning(mfilename, sprintf('Open arrival entry %s has no bound activity; ignored.', ...
            lsn.names{eidx}));
        continue;
    end
    openWiring(end+1) = struct('eidx', eidx, 'firstStep', firstStep, ...
        'exits', [replyExits; terminals]); %#ok<AGROW>
end

%% Pass 2: create classes and reply signals
nsteps = length(stepAidx);
stepClass = cell(nsteps, 1);
stepSignal = cell(nsteps, 1);

for i = 1:nsteps
    refTidx = stepRefTask(i);
    if stepClassOwner(i) ~= i
        % Fork/Join/Router steps carry the job unchanged in the class that entered the fork
        continue;
    end
    if refTidx == 0
        % A step of an open arrival chain travels in an open class.
        stepClass{i} = OpenClass(model, stepName{i});
    elseif stepIsThink(i)
        population = lsn.mult(refTidx);
        stepClass{i} = ClosedClass(model, stepName{i}, population, thinkNode(refTidx));
    else
        stepClass{i} = ClosedClass(model, stepName{i}, 0, thinkNode(refTidx));
    end
end
for i = 1:nsteps
    stepClass{i} = stepClass{stepClassOwner(i)};
end

for i = 1:nsteps
    if stepBlocks(i)
        sig = Signal(model, [stepName{i}, '_Reply'], SignalType.REPLY);
        sig.forJobClass(stepClass{i});
        stepSignal{i} = sig;
    end
end

% Signal.m installs RAND routing at every node; clear it so that link()
% only honours the routes set below.
for i = 1:nsteps
    if ~isempty(stepSignal{i})
        for n = 1:length(model.nodes)
            if ~isa(model.nodes{n}, 'Sink')
                model.nodes{n}.setRouting(stepSignal{i}, RoutingStrategy.DISABLED);
            end
        end
    end
end

% Spawn bindings for phase-2 continuations: each completion of the trigger
% class injects a job of the target class at the same station.
for sp = 1:size(spawnPairs, 1)
    stepClass{spawnPairs(sp, 1)}.setSpawnClass(stepClass{spawnPairs(sp, 2)});
end

% Phase-2 token destructor (closed chains): NEGATIVE signal at a station nothing visits, one per reference task -- see _kb/04-networkstruct.md
ph2DumpNode = [];
ph2DestructorOf = containers.Map('KeyType', 'double', 'ValueType', 'any');
if ~isempty(ph2Exits) && any(ph2Exits(:, 4) > 0)
    ph2DumpNode = Queue(model, 'Ph2Sink', SchedStrategy.FCFS);
    for rft = unique(ph2Exits(ph2Exits(:, 4) > 0, 4).')
        sig = ClosedSignal(model, sprintf('Ph2End_%s', lsn.names{rft}), ...
            SignalType.NEGATIVE, thinkNode(rft));
        for n = 1:length(model.nodes)
            if ~isa(model.nodes{n}, 'Sink')
                model.nodes{n}.setRouting(sig, RoutingStrategy.DISABLED);
            end
        end
        ph2DumpNode.setService(sig, Immediate());
        ph2DestructorOf(rft) = sig;
    end
end

%% Pass 3: service times
for i = 1:nsteps
    if ~isempty(stepNode{i})
        % A Router-hosted merge step owns a class, declared Immediate at the reference think delay (as for signals)
        if stepClassOwner(i) == i && isa(stepNode{i}, 'Router')
            if stepRefTask(i) == 0
                % Open chain: no think delay exists; declare at the caller's host station, which the class never visits
                hostStation{stepHost(i)}.setService(stepClass{i}, Immediate());
            else
                tn_ = thinkNode(stepRefTask(i));
                tn_.setService(stepClass{i}, Immediate());
            end
        end
        continue;
    end
    if stepIsThink(i)
        refTidx = stepRefTask(i);
        thinkDist = lsn.think{refTidx};
        tnode = thinkNode(refTidx);
        if isempty(thinkDist) || isa(thinkDist, 'Immediate') || thinkDist.getMean() < GlobalConstants.FineTol
            tnode.setService(stepClass{i}, Immediate());
        else
            tnode.setService(stepClass{i}, thinkDist);
        end
    else
        station = hostStation{stepHost(i)};
        if isempty(stepSvc{i})
            station.setService(stepClass{i}, Immediate());
        else
            station.setService(stepClass{i}, stepSvc{i});
        end
    end
end

% A reply signal is consumed at the caller's station; declared at every station since a signal with no pending reply falls back to its reference station
for i = 1:nsteps
    if ~isempty(stepSignal{i})
        for h = 1:lsn.nhosts
            hostStation{h}.setService(stepSignal{i}, Immediate());
        end
        tkeys = cell2mat(thinkNode.keys);
        for kk = tkeys
            tn = thinkNode(kk);
            tn.setService(stepSignal{i}, Immediate());
        end
    end
end

%% Cache read/hit/miss wiring, now that the classes exist
for w = 1:length(cacheWiring)
    cw = cacheWiring(w);
    cw.node.setReadItemEntry(stepClass{cw.readStep}, cw.itemproc, cw.nitems);
    cw.node.setHitClass(stepClass{cw.readStep}, stepClass{cw.hitStep});
    cw.node.setMissClass(stepClass{cw.readStep}, stepClass{cw.missStep});
end

%% Pass 4: routing
P = model.initRoutingMatrix();

for e = 1:size(flow, 1)
    i = flow(e, 1); j = flow(e, 2); p = flow(e, 3);
    if flow(e, 5)
        P{stepClass{j}, stepClass{j}}(stationOf(i), stationOf(j)) = p;
    elseif flow(e, 4)
        P{stepSignal{i}, stepClass{j}}(stationOf(i), stationOf(j)) = p;
    else
        P{stepClass{i}, stepClass{j}}(stationOf(i), stationOf(j)) = p;
    end
end

for e = 1:size(reply, 1)
    i = reply(e, 1); owner = reply(e, 2); viaSignal = reply(e, 3); p = reply(e, 4);
    if viaSignal
        % A nested call returns through its own reply signal, which
        % class-switches into the reply signal of the outer call site.
        P{stepSignal{i}, stepSignal{owner}}(stationOf(i), stationOf(owner)) = p;
    else
        P{stepClass{i}, stepSignal{owner}}(stationOf(i), stationOf(owner)) = p;
    end
end

%% Phase-2 chain ends: destroy the spawned token
for e = 1:size(ph2Exits, 1)
    i = ph2Exits(e, 1); viaSig = ph2Exits(e, 2); p = ph2Exits(e, 3); rft = ph2Exits(e, 4);
    if rft == 0
        % Open chain: the spawned token leaves through the Sink.
        ecls = stepClass{i};
        P{ecls, ecls}(stationOf(i), snkNode) = p;
    elseif viaSig
        P{stepSignal{i}, ph2DestructorOf(rft)}(stationOf(i), ph2DumpNode) = p;
    else
        P{stepClass{i}, ph2DestructorOf(rft)}(stationOf(i), ph2DumpNode) = p;
    end
end

%% Open arrival wiring: Source into the first step, exits into the Sink
for w = 1:length(openWiring)
    ow = openWiring(w);
    firstCls = stepClass{ow.firstStep};
    srcNode.setArrival(firstCls, lsn.arrival{ow.eidx});
    P{firstCls, firstCls}(srcNode, stationOf(ow.firstStep)) = 1.0;
    for s = 1:size(ow.exits, 1)
        % Open chains carry no signals, so every exit is an ordinary class.
        ecls = stepClass{ow.exits(s, 1)};
        P{ecls, ecls}(stationOf(ow.exits(s, 1)), snkNode) = ow.exits(s, 3);
    end
end

model.link(P);

%% Thread pools: one finite capacity region, one linear constraint per task
% One admission row A(t,:)*x <= mult(t) per task, coefficients shared across nested tasks -- see _kb/04-networkstruct.md
fcrList = find(fcrTask(:).');
if ~isempty(fcrList)
    Amat = zeros(length(fcrList), length(model.classes));
    bvec = zeros(length(fcrList), 1);
    regionNodes = {};
    for tsel = 1:length(fcrList)
        for stp = 1:nsteps
            if isempty(stepTasks{stp}) || ~any(stepTasks{stp} == fcrList(tsel))
                continue;
            end
            Amat(tsel, stepClass{stp}.index) = 1;
            nd = stationOf(stp);
            if isa(nd, 'Station') && ~any(cellfun(@(x) x == nd, regionNodes))
                regionNodes{end+1} = nd; %#ok<AGROW>
            end
        end
        bvec(tsel) = lsn.mult(fcrList(tsel));
    end
    if any(Amat(:)) && ~isempty(regionNodes)
        fcr = model.addRegion(regionNodes);
        fcr.setConstraint(Amat, bvec);
    end
end

%% ---------------------------------------------------------------- helpers

    function node = stationOf(i)
        if ~isempty(stepNode{i})
            node = stepNode{i};
        elseif stepIsThink(i)
            node = thinkNode(stepRefTask(i));
        else
            node = hostStation{stepHost(i)};
        end
    end

    function id = addStep(aidx, hidx, svc, name, blocks, isthink, refTidx)
        stepAidx(end+1) = ifempty(aidx, 0); %#ok<AGROW>
        stepHost(end+1) = ifempty(hidx, 0); %#ok<AGROW>
        stepSvc{end+1} = svc; %#ok<AGROW>
        stepName{end+1} = name; %#ok<AGROW>
        stepBlocks(end+1) = blocks; %#ok<AGROW>
        stepIsThink(end+1) = isthink; %#ok<AGROW>
        stepRefTask(end+1) = refTidx; %#ok<AGROW>
        stepNode{end+1} = []; %#ok<AGROW>
        id = length(stepAidx);
        stepClassOwner(end+1) = id; %#ok<AGROW>
        % Every thread-pool task on the stack holds a thread here (a synchronous caller releases it only on reply)
        if isempty(threadStack)
            stepTasks{end+1} = []; %#ok<AGROW>
        else
            stepTasks{end+1} = unique(threadStack(fcrTask(threadStack))); %#ok<AGROW>
        end
    end

    function id = addAuxStep(nodeObj, ownerStep, name, refTidx)
        % A step on a Fork, Join or Router node: no station, no service, and
        % no class of its own.
        id = addStep([], [], [], name, false, false, refTidx);
        stepNode{id} = nodeObj;
        stepClassOwner(id) = stepClassOwner(ownerStep);
    end

    function [firstStep, replyExits, terminals] = expandEntry(eidx, refTidx)
        % Expands the activity subgraph bound to an entry, in the call
        % context given by the current entry stack.
        firstStep = [];
        replyExits = zeros(0, 3);
        terminals = zeros(0, 3);

        if any(entryStack == eidx)
            line_warning(mfilename, sprintf('Recursive call cycle at entry %s truncated.', lsn.names{eidx}));
            return;
        end
        entryStack(end+1) = eidx;
        threadStack(end+1) = lsn.parent(eidx);
        restore = onCleanup(@() popEntry());

        if eidx > length(lsn.actsof) || isempty(lsn.actsof{eidx})
            return;
        end

        % The bound activity is the activity successor of the entry.
        succ = find(lsn.graph(eidx, :));
        boundActs = succ(arrayfun(@(a) a <= length(lsn.type) && lsn.type(a) == LayeredNetworkElement.ACTIVITY, succ));
        boundActs = intersect(boundActs, lsn.actsof{eidx}, 'stable');
        if isempty(boundActs)
            return;
        end

        [firstStep, replyExits, terminals] = expandActivities(boundActs(1), eidx, refTidx);

        % Forwarding: with prob p the entry hands off to another entry, which replies directly to the original caller
        fwd = forwardingOf(eidx);
        if ~isempty(fwd) && ~isempty(replyExits)
            ownPorts = replyExits;
            fwdExits = zeros(0, 3);
            pforw = 0.0;
            % Forwarder's thread released at handoff; the forwarded chain expands without it on the thread stack
            fwdThread = threadStack(end);
            threadStack(end) = [];
            for f = 1:size(fwd, 1)
                [fFirst, fReplies, fTerms] = expandEntry(fwd(f, 1), refTidx);
                if isempty(fFirst)
                    continue;
                end
                p = fwd(f, 2);
                for r = 1:size(ownPorts, 1)
                    addRoute(ownPorts(r, 1:2), fFirst, ownPorts(r, 3) * p);
                end
                pforw = pforw + p;
                fwdExits = [fwdExits; fReplies; fTerms]; %#ok<AGROW>
            end
            threadStack(end+1) = fwdThread;
            % What is left of each of this entry's own ports still replies.
            replyExits(:, 3) = replyExits(:, 3) * max(0.0, 1.0 - pforw);
            replyExits = [replyExits; fwdExits];
        end
    end

    function fwd = forwardingOf(eidx)
        % Rows [targetEidx, probability] of the forwarding calls of an entry.
        fwd = zeros(0, 2);
        if ~isfield(lsn, 'calltype') || isempty(lsn.calltype)
            return;
        end
        for cidx = 1:size(lsn.callpair, 1)
            if full(lsn.calltype(cidx)) ~= CallType.FWD || lsn.callpair(cidx, 1) ~= eidx
                continue;
            end
            p = 1.0;
            if isfield(lsn, 'callproc') && ~isempty(lsn.callproc) && ...
                    cidx <= length(lsn.callproc) && isa(lsn.callproc{cidx}, 'Distribution')
                p = lsn.callproc{cidx}.getMean();
            end
            p = min(max(p, 0.0), 1.0);
            if p > GlobalConstants.FineTol
                fwd(end+1, :) = [lsn.callpair(cidx, 2), p]; %#ok<AGROW>
            end
        end
    end

    function popEntry()
        entryStack(end) = [];
        threadStack(end) = [];
    end

    function [firstStep, replyExits, terminals] = expandActivities(a0, eidx, refTidx)
        % Walks the intra-task activity graph from a0, creating steps and
        % expanding every synchronous call site. Reply exits and terminals
        % are returned as ports, [step, isSignal] rows.
        firstStep = [];
        replyExits = zeros(0, 3);
        terminals = zeros(0, 3);

        tidx = lsn.parent(a0);
        localActs = lsn.actsof{eidx};
        visited = containers.Map('KeyType', 'double', 'ValueType', 'any'); % aidx -> [entryStep, exitStep, exitIsSignal]
        joinOf = containers.Map('KeyType', 'double', 'ValueType', 'any');  % join activity -> [joinStep, entryStep]
        forkOwnerStack = [];   % class-owner step of each enclosing AND-fork
        sawReply = false;

        [firstStep, ~] = walk(a0);

        % Every chain end returns the token to the caller: as a reply if the entry replies anywhere, else a plain terminal
        if sawReply
            replyExits = [replyExits; terminals]; %#ok<AGROW>
            terminals = zeros(0, 3);
        end

        function [entryStep, exitPort] = walk(aidx)
            if isKey(visited, aidx)
                se = visited(aidx);
                entryStep = se(1);
                exitPort = se(2:3);
                return;
            end
            % The activity bound to an ItemEntry of a CacheTask is the read
            % step: it sits on the Cache node rather than on the processor.
            cacheNode = [];
            if tidx <= length(lsn.iscache) && lsn.iscache(tidx) && full(lsn.graph(eidx, aidx)) > 0
                cacheNode = getCacheNode(tidx);
            end
            [entryStep, exitPort] = makeActivitySteps(aidx, tidx, refTidx, cacheNode);
            visited(aidx) = [entryStep, exitPort];

            % Reply is deferred to the ends of the chain, not emitted here -- see _kb/04-networkstruct.md
            repliesHere = false;
            if isfield(lsn, 'replygraph') && ~isempty(lsn.replygraph)
                a = aidx - lsn.ashift;
                e = eidx - lsn.eshift;
                if a >= 1 && a <= size(lsn.replygraph, 1) && e >= 1 && e <= size(lsn.replygraph, 2)
                    repliesHere = full(lsn.replygraph(a, e));
                end
            end
            if repliesHere
                sawReply = true;
            end

            % Local successors within the same entry.
            succ = find(lsn.graph(aidx, :));
            succ = succ(ismember(succ, localActs));
            if isempty(succ)
                terminals(end+1, :) = [exitPort, 1.0]; %#ok<AGROW>
                return;
            end

            % Phase 2 runs after the reply via sn.classspawn; needs a station-departure exit in the step's own class -- see _kb/04-networkstruct.md
            if repliesHere
                if ~isempty(cacheNode) && length(succ) >= 2
                    % Phase 2 at a cache read: each hit/miss outcome routes through its own immediate trigger step -- see _kb/04-networkstruct.md
                    trigH = addStep(aidx, lsn.parent(tidx), [], ...
                        [lsn.names{aidx}, '_ph2h'], false, false, refTidx);
                    trigM = addStep(aidx, lsn.parent(tidx), [], ...
                        [lsn.names{aidx}, '_ph2m'], false, false, refTidx);
                    addCacheRoute(entryStep, trigH);
                    addCacheRoute(entryStep, trigM);
                    cacheWiring(end+1) = struct('node', cacheNode, ...
                        'readStep', entryStep, 'hitStep', trigH, 'missStep', trigM, ...
                        'itemproc', lsn.itemproc{eidx}, 'nitems', lsn.nitems(eidx)); %#ok<AGROW>
                    replyExits(end+1, :) = [trigH, 0, 1.0]; %#ok<AGROW>
                    replyExits(end+1, :) = [trigM, 0, 1.0]; %#ok<AGROW>
                    savedStack = threadStack;
                    threadStack = tidx;
                    nT0 = size(terminals, 1);
                    for hm = 1:2
                        [sEntry, ~] = walk(succ(hm));
                        if ~isempty(stepNode{sEntry})
                            hmHead = addStep(aidx, lsn.parent(tidx), [], ...
                                sprintf('%s_ph2b%d', lsn.names{aidx}, hm), false, false, refTidx);
                            addRoute([hmHead, 0], sEntry, 1.0);
                            sEntry = hmHead;
                        end
                        if hm == 1
                            spawnPairs(end+1, :) = [trigH, sEntry]; %#ok<AGROW>
                        else
                            spawnPairs(end+1, :) = [trigM, sEntry]; %#ok<AGROW>
                        end
                    end
                    ph2New = terminals(nT0+1:end, :);
                    terminals(nT0+1:end, :) = [];
                    for r2 = 1:size(ph2New, 1)
                        ph2Exits(end+1, :) = [ph2New(r2, 1:3), refTidx]; %#ok<AGROW>
                    end
                    threadStack = savedStack;
                    return;
                end
                % Inside an AND-fork branch, the lift applies only at the branch tail -- see _kb/04-networkstruct.md
                okCtx = isempty(cacheNode) && ...
                    (isempty(forkOwnerStack) || isAndJoinPre(aidx));
                % A merge step at the host station normalises a call-site exit to a station departure before the phase-2 walk
                if okCtx && (exitPort(2) == 1 || ...
                        (~isempty(stepNode{exitPort(1)}) && isa(stepNode{exitPort(1)}, 'Router')))
                    trig = addStep(aidx, lsn.parent(tidx), [], ...
                        [lsn.names{aidx}, '_ph2t'], false, false, refTidx);
                    addRoute(exitPort, trig, 1.0);
                    exitPort = [trig, 0];
                end
                ph2Spawn = okCtx && exitPort(2) == 0 && isempty(stepNode{exitPort(1)});
                if ~ph2Spawn
                    line_warning(mfilename, sprintf(['Phase-2 activities of %s run ' ...
                        'before the reply: the boundary is not a station ' ...
                        'departure, a degenerate cache read, or mid-branch ' ...
                        'inside an AND-fork.'], lsn.names{aidx}));
                else
                    replyExits(end+1, :) = [exitPort, 1.0]; %#ok<AGROW>
                    savedStack = threadStack;
                    threadStack = tidx;
                    nT0 = size(terminals, 1);
                    posSucc = succ(full(lsn.graph(aidx, succ)) > 0);
                    target = [];
                    if isAndFork(succ)
                        % Phase 2 opens with an AND-fork: spawn into an immediate head step and fork from there
                        head = addStep(aidx, lsn.parent(tidx), [], ...
                            [lsn.names{aidx}, '_ph2'], false, false, refTidx);
                        wireAndFork([head, 0], succ, aidx);
                        target = head;
                    elseif isAndJoinPre(aidx)
                        % Phase 2 at an AND-join branch tail: spawn into an immediate head standing in for this branch at the Join
                        head = addStep(aidx, lsn.parent(tidx), [], ...
                            [lsn.names{aidx}, '_ph2'], false, false, refTidx);
                        wireAndJoin([head, 0], succ(1));
                        target = head;
                    elseif isscalar(posSucc)
                        [sEntry, ~] = walk(posSucc);
                        if isempty(stepNode{sEntry})
                            target = sEntry;
                        end
                    end
                    if isempty(target)
                        % Branching phase 2, or a head on a non-station node: spawn into an immediate head carrying branch probabilities
                        head = addStep(aidx, lsn.parent(tidx), [], ...
                            [lsn.names{aidx}, '_ph2'], false, false, refTidx);
                        for s2 = posSucc
                            [sEntry, ~] = walk(s2);
                            addRoute([head, 0], sEntry, full(lsn.graph(aidx, s2)));
                        end
                        target = head;
                    end
                    spawnPairs(end+1, :) = [exitPort(1), target]; %#ok<AGROW>
                    ph2New = terminals(nT0+1:end, :);
                    terminals(nT0+1:end, :) = [];
                    for r2 = 1:size(ph2New, 1)
                        ph2Exits(end+1, :) = [ph2New(r2, 1:3), refTidx]; %#ok<AGROW>
                    end
                    threadStack = savedStack;
                    return;
                end
            end

            if ~isempty(cacheNode)
                % CacheAccess precedence: successors are hit then miss branch; the Cache node decides, so no probability on these routes
                if length(succ) < 2
                    line_warning(mfilename, sprintf(...
                        'Cache read %s has no hit/miss pair; treated as an ordinary activity.', ...
                        lsn.names{aidx}));
                else
                    [hEntry, ~] = walk(succ(1));
                    [mEntry, ~] = walk(succ(2));
                    addCacheRoute(entryStep, hEntry);
                    addCacheRoute(entryStep, mEntry);
                    cacheWiring(end+1) = struct('node', cacheNode, ...
                        'readStep', entryStep, 'hitStep', hEntry, 'missStep', mEntry, ...
                        'itemproc', lsn.itemproc{eidx}, 'nitems', lsn.nitems(eidx)); %#ok<AGROW>
                    return;
                end
            end

            if isAndFork(succ)
                wireAndFork(exitPort, succ, aidx);
                return;
            end

            if isAndJoinPre(aidx)
                wireAndJoin(exitPort, succ(1));
                return;
            end

            for s = succ
                p = full(lsn.graph(aidx, s));
                if p <= 0
                    continue;
                end
                [sEntry, ~] = walk(s);
                addRoute(exitPort, sEntry, p);
            end
        end

        function wireAndFork(fromPort, fsucc, aidx)
            % AND-fork: the branches run concurrently, so the job is
            % replicated by a Fork node. One Router per branch, because a
            % Fork cannot switch class per output link, and the branch
            % class is what tells the branch apart.
            forkNode = Fork(model, ['Fork_', lsn.names{aidx}]);
            forkStep = addAuxStep(forkNode, fromPort(1), ['Fork_', lsn.names{aidx}], refTidx);
            addRoute(fromPort, forkStep, 1.0);
            forkOwnerStack(end+1) = forkStep;
            % Walk a replying branch first so the Join/post-join subgraph is created in its phase-2 context
            rep = arrayfun(@branchReplies, fsucc);
            fsucc = [fsucc(rep), fsucc(~rep)];
            for b = 1:length(fsucc)
                routerNode = Router(model, sprintf('Fork_%s_%d', lsn.names{aidx}, b));
                routerStep = addAuxStep(routerNode, forkStep, ...
                    sprintf('Fork_%s_%d', lsn.names{aidx}, b), refTidx);
                addRoute([forkStep, 0], routerStep, 1.0);
                [sEntry, ~] = walk(fsucc(b));
                addRoute([routerStep, 0], sEntry, 1.0);
            end
            forkOwnerStack(end) = [];
        end

        function wireAndJoin(fromPort, joinAidx)
            % Route a branch tail into the AND-join, creating the Join on
            % first arrival. The siblings are matched there in the class
            % that entered the fork, so switch back to it on the way in.
            % The post-join subgraph is walked in the calling context.
            if isKey(joinOf, joinAidx)
                js = joinOf(joinAidx);
                addRoute(fromPort, js(1), 1.0);
                return;
            end
            if isempty(forkOwnerStack)
                line_warning(mfilename, sprintf(...
                    'AND-join at %s has no enclosing AND-fork; branches are serialised.', ...
                    lsn.names{joinAidx}));
                [sEntry, ~] = walk(joinAidx);
                addRoute(fromPort, sEntry, 1.0);
                return;
            end
            joinNode = Join(model, ['Join_', lsn.names{joinAidx}], stepNode{forkOwnerStack(end)});
            joinStep = addAuxStep(joinNode, forkOwnerStack(end), ...
                ['Join_', lsn.names{joinAidx}], refTidx);
            addRoute(fromPort, joinStep, 1.0);
            [sEntry, ~] = walk(joinAidx);
            addRoute([joinStep, 0], sEntry, 1.0);
            joinOf(joinAidx) = [joinStep, sEntry];
        end

        function tf = branchReplies(a0)
            % True if the branch rooted at a0 contains an activity that
            % replies to the current entry, searching up to and excluding
            % the AND-join that closes the branch.
            tf = false;
            stack = a0;
            seenActs = [];
            while ~isempty(stack)
                a = stack(end);
                stack(end) = [];
                if any(seenActs == a)
                    continue;
                end
                seenActs(end+1) = a; %#ok<AGROW>
                if isfield(lsn, 'replygraph') && ~isempty(lsn.replygraph)
                    ar = a - lsn.ashift;
                    er = eidx - lsn.eshift;
                    if ar >= 1 && ar <= size(lsn.replygraph, 1) && ...
                            er >= 1 && er <= size(lsn.replygraph, 2) && ...
                            full(lsn.replygraph(ar, er))
                        tf = true;
                        return;
                    end
                end
                if isAndJoinPre(a)
                    continue;
                end
                nxt = find(lsn.graph(a, :));
                stack = [stack, nxt(ismember(nxt, localActs))]; %#ok<AGROW>
            end
        end
    end

    function [entryStep, exitPort] = makeActivitySteps(aidx, tidx, refTidx, cacheNode)
        % One step for the host demand, plus one step per unrolled
        % synchronous call stage.
        hidx = lsn.parent(tidx);
        if nargin >= 4 && ~isempty(cacheNode)
            % A read step holds no demand and issues no call: the lookup is
            % instantaneous and the work is done on the hit or miss branch.
            entryStep = addStep(aidx, hidx, [], lsn.names{aidx}, false, false, refTidx);
            stepNode{entryStep} = cacheNode;
            if aidx <= length(lsn.callsof) && ~isempty(lsn.callsof{aidx})
                line_warning(mfilename, sprintf(...
                    'Calls issued by cache read activity %s are ignored.', lsn.names{aidx}));
            end
            if aidx <= length(lsn.hostdem) && isa(lsn.hostdem{aidx}, 'Distribution') && ...
                    ~isa(lsn.hostdem{aidx}, 'Immediate') && lsn.hostdem{aidx}.getMean() > GlobalConstants.FineTol
                line_warning(mfilename, sprintf(...
                    'Host demand of cache read activity %s is ignored.', lsn.names{aidx}));
            end
            exitPort = [entryStep, 0];
            return;
        end
        svc = [];
        if aidx <= length(lsn.hostdem) && isa(lsn.hostdem{aidx}, 'Distribution')
            d = lsn.hostdem{aidx};
            if ~isa(d, 'Immediate') && d.getMean() > GlobalConstants.FineTol
                svc = d;
            end
        end

        callStages = synchCallStages(aidx);
        entryStep = addStep(aidx, hidx, svc, lsn.names{aidx}, false, false, refTidx);
        % cur is the port through which the activity is currently left. A
        % blocking call site is left through its reply signal.
        cur = [entryStep, 0];

        % A call blocks the caller's server only when host servers are finite, the task is not a thread pool,
        % multiplicity is finite, and the chain is not open -- see _kb/04-networkstruct.md
        hostBlocks = ~hostIsDelay(hidx) && ~fcrTask(tidx) && refTidx ~= 0 && ...
            isfinite(lsn.mult(tidx)) && lsn.sched(tidx) ~= SchedStrategy.INF;

        for k = 1:length(callStages)
            stage = callStages(k);
            % Async call does not hold the caller's server (returns via an external class switch); still serialised, which is the approximation
            blocks = hostBlocks && ~stage.isasync;
            % Async send also releases the caller's thread for the callee expansion
            if stage.isasync
                asyncThread = threadStack(end);
                threadStack(end) = [];
            end
            [calleeFirst, calleeReplies, calleeTerms] = expandEntry(stage.targetEidx, refTidx);
            if stage.isasync
                threadStack(end+1) = asyncThread;
            end
            if isempty(calleeFirst)
                continue;   % callee not expandable: drop the call, never block
            end
            % A callee path that neither replies nor continues still holds a
            % token, so it returns to the caller like a reply would.
            calleeReplies = [calleeReplies; calleeTerms]; %#ok<AGROW>

            % Merge step needed when the call may be skipped, another stage follows, or an AND-join branch tail must reach
            % the Join in an ordinary class (a REPLY signal carries no forked-task identity) -- see _kb/04-networkstruct.md
            needsMerge = (stage.prob < 1.0) || (k < length(callStages)) || ...
                ~blocks || isAndJoinPre(aidx);
            if needsMerge
                nxt = addStep(aidx, hidx, [], sprintf('%s_c%d_ret', lsn.names{aidx}, k), false, false, refTidx);
                if ~blocks
                    % Non-blocking return carries no signal, so the merge point can sit on a Router (no server to queue behind)
                    stepNode{nxt} = Router(model, sprintf('%s_c%d_ret', lsn.names{aidx}, k));
                end
            end

            if blocks
                % First mandatory call bound to the service class itself; an intervening class switch would release the server
                if k == 1 && stage.prob >= 1.0 && isequal(cur, [entryStep, 0])
                    blk = entryStep;
                    stepBlocks(blk) = true;
                else
                    blk = addStep(aidx, hidx, [], sprintf('%s_c%d', lsn.names{aidx}, k), true, false, refTidx);
                    addRoute(cur, blk, stage.prob);
                    if stage.prob < 1.0
                        addRoute(cur, nxt, 1.0 - stage.prob);
                    end
                end
                addRoute([blk, 0], calleeFirst, 1.0);
                for r = 1:size(calleeReplies, 1)
                    reply(end+1, :) = [calleeReplies(r, 1), blk, calleeReplies(r, 2), calleeReplies(r, 3)]; %#ok<AGROW>
                end
                if needsMerge
                    addRoute([blk, 1], nxt, 1.0);
                    cur = [nxt, 0];
                else
                    cur = [blk, 1];
                end
            else
                addRoute(cur, calleeFirst, stage.prob);
                if stage.prob < 1.0
                    addRoute(cur, nxt, 1.0 - stage.prob);
                end
                for r = 1:size(calleeReplies, 1)
                    addRoute(calleeReplies(r, 1:2), nxt, calleeReplies(r, 3));
                end
                cur = [nxt, 0];
            end
        end
        exitPort = cur;
    end

    function cnode = getCacheNode(tidx)
        % One Cache node per CacheTask. The name is suffixed so that it does
        % not collide with the station of the task's processor.
        if isKey(cacheNodeOf, tidx)
            cnode = cacheNodeOf(tidx);
            return;
        end
        cnode = Cache(model, [lsn.names{tidx}, '_Cache'], lsn.nitems(tidx), ...
            lsn.itemcap{tidx}, lsn.replacestrat(tidx));
        cacheNodeOf(tidx) = cnode;
    end

    function tf = isAndFork(succ)
        % An AND-fork is a precedence whose post activities are all marked
        % POST_AND: they are entered together, not with a probability each.
        tf = false;
        if length(succ) < 2 || ~isfield(lsn, 'actposttype') || isempty(lsn.actposttype)
            return;
        end
        tf = all(arrayfun(@(s) s <= length(lsn.actposttype) && ...
            full(lsn.actposttype(s)) == ActivityPrecedenceType.POST_AND, succ));
    end

    function tf = isAndJoinPre(aidx)
        % An activity marked PRE_AND is one branch tail of an AND-join.
        tf = isfield(lsn, 'actpretype') && ~isempty(lsn.actpretype) && ...
            aidx <= length(lsn.actpretype) && ...
            full(lsn.actpretype(aidx)) == ActivityPrecedenceType.PRE_AND;
    end

    function addRoute(fromPort, toStep, prob)
        flow(end+1, :) = [fromPort(1), toStep, prob, fromPort(2), 0]; %#ok<AGROW>
    end

    function addCacheRoute(fromStep, toStep)
        % Leaving a Cache node: the switch into the hit or the miss class is
        % made by the node, so the route is declared in the target class.
        flow(end+1, :) = [fromStep, toStep, 1.0, 0, 1]; %#ok<AGROW>
    end

    function stages = synchCallStages(aidx)
        % Unrolls the synchronous calls of an activity into blocking stages.
        stages = struct('targetEidx', {}, 'prob', {}, 'isasync', {});
        if aidx > length(lsn.callsof) || isempty(lsn.callsof{aidx})
            return;
        end
        for cidx = lsn.callsof{aidx}
            ctype = full(lsn.calltype(cidx));
            if ctype ~= CallType.SYNC && ctype ~= CallType.ASYNC
                continue;
            end
            isasync = (ctype == CallType.ASYNC);
            targetEidx = lsn.callpair(cidx, 2);
            m = 1.0;
            if isfield(lsn, 'callproc') && ~isempty(lsn.callproc) && ...
                    cidx <= length(lsn.callproc) && isa(lsn.callproc{cidx}, 'Distribution')
                m = lsn.callproc{cidx}.getMean();
            end
            nfull = floor(m + GlobalConstants.FineTol);
            frac = m - nfull;
            if nfull > MAXCALLSTAGES
                line_warning(mfilename, sprintf('Call multiplicity %g on %s truncated to %d stages.', ...
                    m, lsn.callnames{cidx}, MAXCALLSTAGES));
                nfull = MAXCALLSTAGES;
                frac = 0;
            end
            for k = 1:nfull
                stages(end+1) = struct('targetEidx', targetEidx, 'prob', 1.0, 'isasync', isasync); %#ok<AGROW>
            end
            if frac > GlobalConstants.FineTol
                stages(end+1) = struct('targetEidx', targetEidx, 'prob', frac, 'isasync', isasync); %#ok<AGROW>
            end
        end
    end

    function warnUnsupported(lsn)
        if any(full(lsn.calltype) == CallType.ASYNC)
            line_warning(mfilename, ['Asynchronous calls are represented by LQN2QN as ' ...
                'non-blocking visits: the caller releases its server but remains ' ...
                'serialised behind the callee.']);
        end
        if isfield(lsn, 'hasretrieval') && any(lsn.hasretrieval)
            line_warning(mfilename, 'Delayed-hit retrieval on the cache miss path is not represented by LQN2QN.');
        end
    end

    function tf = taskHasAndFork(tidx_)
        % True if any activity of the task is the head of an AND-fork
        % branch, i.e. is marked POST_AND.
        tf = false;
        if ~isfield(lsn, 'actposttype') || isempty(lsn.actposttype)
            return;
        end
        for a_ = lsn.ashift+1:lsn.ashift+lsn.nacts
            if lsn.parent(a_) == tidx_ && a_ <= length(lsn.actposttype) && ...
                    full(lsn.actposttype(a_)) == ActivityPrecedenceType.POST_AND
                tf = true;
                return;
            end
        end
    end

    function v = ifempty(x, d)
        if isempty(x)
            v = d;
        else
            v = x;
        end
    end

end
