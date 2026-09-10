function model = LQN2QN(lqn, replication)
% LQN2QN Convert a LayeredNetwork (LQN) to a Network (QN) using REPLY signals
%
% model = LQN2QN(lqn) flattens a LayeredNetwork into a single queueing
% network in which synchronous call blocking is represented by REPLY
% signals.
%
% model = LQN2QN(lqn, replication) selects how task and processor
% replication is represented: 'auto' (default) materialises the replicas
% while the expansion stays within the instantiation budget and pools them
% otherwise, 'materialize' always materialises, 'pool' always pools.
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
% - An AND-join quorum k of n is applied to the Join node in the class that
%   entered the Fork; k equal to the branch count is the default wait-for-all
%   and is left alone. An activity think time becomes an extra step on a shared
%   ActivityThink delay, in series with the host demand, so the task keeps its
%   thread for it while its processor is released.
%
% - Replication is represented in one of two ways. Under materialisation each
%   replica of a processor is a station of its own and each replica of a task
%   carries its own copy of the expanded step graph, its own reply signals and
%   its own admission row; a call from replica i of the caller reaches the
%   fan-out block {(i*f+k) mod r} of the callee replicas and splits its call
%   mean uniformly over them, which is deterministic pairing at f=1 and a
%   uniform broadcast at f=r. Under pooling the replicas of a processor
%   collapse into one station of r times the servers, a replicated thread pool
%   into one admission row of r times the bound, and a replicated reference
%   task into one class of r times the population; that is exact at an
%   infinite-server host and optimistic elsewhere, since pooled servers share
%   one queue while the replicas hold r separate ones.
%
% - A CacheTask with delayed-hit retrieval gets a retrieval system, which is an
%   ordinary queueing network: one PS fetch station per cache replica, entered
%   and left by the read class, with the Cache node coalescing concurrent misses
%   of the same item. The fetch is what the miss branch does, so the miss
%   activity's host demand moves onto that station; calls issued by the miss
%   activity stay outside the retrieval system and are warned.
%
% - A SetupTask carries its setup and delay-off times onto its host station
%   as the Queue setup/delay-off pair, per step class: the server shuts down
%   after the delay-off idle period and pays the setup on the next arrival. An
%   infinite-server processor never shuts down, so the pair is dropped there
%   with a warning, as is a setup with no delay-off time.
%
% Not yet represented: retrieval on a cache read with phase-2 successors and the
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

if nargin < 2 || isempty(replication)
    replication = 'auto';
end
if ~any(strcmp(replication, {'auto', 'materialize', 'pool'}))
    line_error(mfilename, 'replication must be ''auto'', ''materialize'' or ''pool''.');
end

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
% Above this many replica subgraph instantiations 'auto' pools instead of
% materialising: the routing matrix is dense in (classes x nodes), so the
% conversion cost grows quadratically in the instantiation count.
MAXREPLINSTANCES = 128;

%% Unsupported features, reported once each
warnUnsupported(lsn);

%% Replication
% A replicated element is r identical copies of itself, either materialised one
% station and one step-graph copy per replica, or pooled -- see _kb/06-solver-catalog.md
replRaw = ones(lsn.nhosts + lsn.ntasks, 1);
if isfield(lsn, 'repl') && ~isempty(lsn.repl)
    raw_ = full(lsn.repl(:));
    nr_ = min(length(raw_), lsn.nhosts + lsn.ntasks);
    replRaw(1:nr_) = max(1, round(raw_(1:nr_)));
end
replicated = find(replRaw > 1).';
% Stride of the (element, replica) composite keys used for maps and stacks
RKEY = max(1, max(replRaw));

materialize = ~isempty(replicated) && (strcmp(replication, 'materialize') || ...
    (strcmp(replication, 'auto') && replInstantiations() <= MAXREPLINSTANCES));
if ~isempty(replicated) && ~materialize
    line_warning(mfilename, sprintf(['Replication of %s is pooled: its replicas become ' ...
        'one station of r times the servers, one admission row of r times the bound ' ...
        'and one reference class of r times the population, which is exact at an ' ...
        'infinite-server host and optimistic elsewhere. Pass ''materialize'' for one ' ...
        'station and one step-graph copy per replica.'], lsn.names{replicated(1)}));
end

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

%% Stations: one per host processor replica
hostStation = cell(lsn.nhosts, RKEY);
hostIsDelay = false(lsn.nhosts, RKEY);
for h = 1:lsn.nhosts
    nservers = lsn.mult(h);
    sched = lsn.sched(h);
    for m_ = 0:nrep(h)-1
        if isinf(nservers) || sched == SchedStrategy.INF
            hostStation{h, m_+1} = Delay(model, suffixed(lsn.names{h}, m_));
            hostIsDelay(h, m_+1) = true;
        else
            q = Queue(model, suffixed(lsn.names{h}, m_), sched);
            % A pooled processor carries the servers of all its replicas
            q.setNumberOfServers(nservers * poolFactor(h));
            hostStation{h, m_+1} = q;
        end
    end
end

%% Think delays, one per reference task replica
thinkNode = configureDictionary('double','cell');
for rt = 1:length(refTaskIndices)
    refTidx = refTaskIndices(rt);
    for m_ = 0:nrep(refTidx)-1
        thinkNode{ekey(refTidx, m_)} = Delay(model, ...
            suffixed([lsn.names{refTidx}, '_Think'], m_));
    end
end

%% Pass 1: expand the activity graph into a step graph
% Steps carry no LINE objects yet, so all classes can be created before the routing matrix is initialised
stepAidx = [];        % activity index of the step
stepHost = zeros(0, 2); % [host processor, replica] of the step station
stepSvc = {};         % service distribution at the step station, [] if none
stepName = {};        % class name
stepBlocks = [];      % true if the step blocks on a synchronous call
stepIsThink = [];     % true for a think step (reference task)
stepRefTask = zeros(0, 2); % [reference task, replica] the step belongs to, [0 0] on an open chain
stepNode = {};        % Fork/Join/Router node the step sits on, [] for stations
stepClassOwner = [];  % step whose class this step travels in (itself, normally)
stepTasks = {};       % thread-pool task replicas holding a thread while at this step, as composite keys

% flow/reply/spawnPairs/ph2Exits step-graph arrays -- see _kb/04-networkstruct.md
flow = zeros(0, 5);
reply = zeros(0, 4);
spawnPairs = zeros(0, 2);
ph2Exits = zeros(0, 5);

entryStack = zeros(0, 2);   % [entry, replica]; guards against recursive call cycles
threadStack = [];  % task replicas holding a thread during expansion, as composite keys; tracks entryStack except across forwarding (releases the forwarder's thread)

% Cache nodes, one per CacheTask, and the read/hit/miss wiring to apply once
% the classes exist.
cacheNodeOf = configureDictionary('double','cell');
fetchNodeOf = configureDictionary('double','cell');
cacheWiring = struct('node', {}, 'readStep', {}, 'hitStep', {}, 'missStep', {}, ...
    'itemproc', {}, 'nitems', {}, 'fetch', {}, 'fetchSvc', {});

usedNames = configureDictionary('string', 'double');  % class and node names must be unique
joinQuorum = zeros(0, 2);   % [joinStep, joinAidx]; applied once the fork class exists
actThinkNode = [];          % shared INF station carrying the activity think times

% One closed chain per reference task replica: the replicas are separate
% populations that meet only where they share a station.
for rt = 1:length(refTaskIndices)
    refTidx = refTaskIndices(rt);
    for rep = 0:nrep(refTidx)-1
        refKey = [refTidx, rep];
        thinkStep = addStep([], [], [], ...
            uniqueName(suffixed([lsn.names{refTidx}, '_Think'], rep)), false, true, refKey);

        entries = lsn.entriesof{refTidx};
        for eidx = entries
            [firstStep, replyExits, terminals] = expandEntry(eidx, refKey, rep);
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
    % Each replica of the entry's task receives its own arrival stream.
    for rep = 0:nrep(lsn.parent(eidx))-1
        [firstStep, replyExits, terminals] = expandEntry(eidx, [0, 0], rep);
        if isempty(firstStep)
            line_warning(mfilename, sprintf('Open arrival entry %s has no bound activity; ignored.', ...
                lsn.names{eidx}));
            continue;
        end
        openWiring(end+1) = struct('eidx', eidx, 'firstStep', firstStep, ...
            'exits', [replyExits; terminals]); %#ok<AGROW>
    end
end

%% Pass 2: create classes and reply signals
nsteps = length(stepAidx);
stepClass = cell(nsteps, 1);
stepSignal = cell(nsteps, 1);

for i = 1:nsteps
    refTidx = stepRefTask(i, 1);
    if stepClassOwner(i) ~= i
        % Fork/Join/Router steps carry the job unchanged in the class that entered the fork
        continue;
    end
    if refTidx == 0
        % A step of an open arrival chain travels in an open class.
        stepClass{i} = OpenClass(model, stepName{i});
    elseif stepIsThink(i)
        % A pooled reference task holds the population of all its replicas
        population = lsn.mult(refTidx) * poolFactor(refTidx);
        stepClass{i} = ClosedClass(model, stepName{i}, population, thinkNode{ekey(refTidx, stepRefTask(i, 2))});
    else
        stepClass{i} = ClosedClass(model, stepName{i}, 0, thinkNode{ekey(refTidx, stepRefTask(i, 2))});
    end
end
for i = 1:nsteps
    stepClass{i} = stepClass{stepClassOwner(i)};
end

% AND-join quorum, in the class the siblings are matched in
for jq_ = 1:size(joinQuorum, 1)
    applyJoinQuorum(stepNode{joinQuorum(jq_, 1)}, stepClass{joinQuorum(jq_, 1)}, joinQuorum(jq_, 2));
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
ph2DestructorOf = configureDictionary('double','cell');
if ~isempty(ph2Exits) && any(ph2Exits(:, 4) > 0)
    ph2DumpNode = Queue(model, 'Ph2Sink', SchedStrategy.FCFS);
    ph2Ref = unique(ph2Exits(ph2Exits(:, 4) > 0, 4:5), 'rows');
    for r_ = 1:size(ph2Ref, 1)
        rft = ph2Ref(r_, 1); rrep = ph2Ref(r_, 2);
        sig = ClosedSignal(model, suffixed(sprintf('Ph2End_%s', lsn.names{rft}), rrep), ...
            SignalType.NEGATIVE, thinkNode{ekey(rft, rrep)});
        for n = 1:length(model.nodes)
            if ~isa(model.nodes{n}, 'Sink')
                model.nodes{n}.setRouting(sig, RoutingStrategy.DISABLED);
            end
        end
        ph2DumpNode.setService(sig, Immediate());
        ph2DestructorOf{ekey(rft, rrep)} = sig;
    end
end

%% Pass 3: service times
for i = 1:nsteps
    if ~isempty(stepNode{i})
        if stepClassOwner(i) == i && ~isempty(actThinkNode) && stepNode{i} == actThinkNode
            actThinkNode.setService(stepClass{i}, stepSvc{i});
            continue;
        end
        % A Router-hosted merge step owns a class, declared Immediate at the reference think delay (as for signals)
        if stepClassOwner(i) == i && isa(stepNode{i}, 'Router')
            if stepRefTask(i, 1) == 0
                % Open chain: no think delay exists; declare at the caller's host station, which the class never visits
                hostStation{stepHost(i, 1), stepHost(i, 2)+1}.setService(stepClass{i}, Immediate());
            else
                tn_ = thinkNode{ekey(stepRefTask(i, 1), stepRefTask(i, 2))};
                tn_.setService(stepClass{i}, Immediate());
            end
        end
        continue;
    end
    if stepIsThink(i)
        refTidx = stepRefTask(i, 1);
        thinkDist = lsn.think{refTidx};
        tnode = thinkNode{ekey(refTidx, stepRefTask(i, 2))};
        if isempty(thinkDist) || isa(thinkDist, 'Immediate') || thinkDist.getMean() < GlobalConstants.FineTol
            tnode.setService(stepClass{i}, Immediate());
        else
            tnode.setService(stepClass{i}, thinkDist);
        end
    else
        station = hostStation{stepHost(i, 1), stepHost(i, 2)+1};
        if isempty(stepSvc{i})
            station.setService(stepClass{i}, Immediate());
        else
            station.setService(stepClass{i}, stepSvc{i});
        end
    end
end

% A SetupTask is a server that shuts down when idle and pays a setup on the
% next arrival, which is the Queue setup/delay-off pair at its host station
warnedSetup = [];
for i = 1:nsteps
    if ~isempty(stepNode{i}) || stepIsThink(i) || stepAidx(i) == 0
        continue;
    end
    tidx_ = lsn.parent(stepAidx(i));
    [su_, do_, status_] = functionTimesOf(tidx_);
    if strcmp(status_, 'none')
        continue;
    end
    if strcmp(status_, 'nodelayoff')
        if ~any(warnedSetup == tidx_)
            line_warning(mfilename, sprintf(['Setup of setup task %s is not represented: ' ...
                'it has no delay-off time, so its server never shuts down and never sets ' ...
                'up again.'], lsn.names{tidx_}));
            warnedSetup(end+1) = tidx_; %#ok<AGROW>
        end
        continue;
    end
    % Delay subclasses Queue, so the infinite server is tested by hostIsDelay
    if hostIsDelay(stepHost(i, 1), stepHost(i, 2)+1)
        if ~any(warnedSetup == tidx_)
            line_warning(mfilename, sprintf(['Setup of setup task %s is not represented: ' ...
                'its processor is an infinite server, which never shuts down.'], lsn.names{tidx_}));
            warnedSetup(end+1) = tidx_; %#ok<AGROW>
        end
        continue;
    end
    hostStation{stepHost(i, 1), stepHost(i, 2)+1}.setDelayOff(stepClass{i}, su_, do_);
end

% A reply signal is consumed at the caller's station; declared at every station since one with no pending reply falls back to its reference station
for i = 1:nsteps
    if ~isempty(stepSignal{i})
        for h = 1:lsn.nhosts
            for m_ = 0:nrep(h)-1
                hostStation{h, m_+1}.setService(stepSignal{i}, Immediate());
            end
        end
        tkeys = keys(thinkNode);
        for kk = tkeys
            tn = thinkNode{kk};
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
    if ~isempty(cw.fetch)
        % Service and routing of the retrieval system are read off the read class
        if isempty(cw.fetchSvc)
            cw.fetch.setService(stepClass{cw.readStep}, Immediate());
        else
            cw.fetch.setService(stepClass{cw.readStep}, cw.fetchSvc);
        end
        cw.node.setRetrievalSystem(stepClass{cw.readStep}, stepClass{cw.missStep}, cw.fetch);
    end
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

%% Retrieval systems: the read class circulates cache -> fetch -> cache
for w = 1:length(cacheWiring)
    cw = cacheWiring(w);
    if ~isempty(cw.fetch)
        rcls = stepClass{cw.readStep};
        P{rcls, rcls}(cw.node, cw.fetch) = 1.0;
        P{rcls, rcls}(cw.fetch, cw.node) = 1.0;
    end
end

%% Phase-2 chain ends: destroy the spawned token
for e = 1:size(ph2Exits, 1)
    i = ph2Exits(e, 1); viaSig = ph2Exits(e, 2); p = ph2Exits(e, 3);
    rft = ph2Exits(e, 4); rrep = ph2Exits(e, 5);
    if rft == 0
        % Open chain: the spawned token leaves through the Sink.
        ecls = stepClass{i};
        P{ecls, ecls}(stationOf(i), snkNode) = p;
    elseif viaSig
        P{stepSignal{i}, ph2DestructorOf{ekey(rft, rrep)}}(stationOf(i), ph2DumpNode) = p;
    else
        P{stepClass{i}, ph2DestructorOf{ekey(rft, rrep)}}(stationOf(i), ph2DumpNode) = p;
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

%% Thread pools: one finite capacity region, one linear constraint per task replica
% One admission row A(t,:)*x <= mult(t) per task replica, coefficients shared across nested tasks -- see _kb/04-networkstruct.md
fcrList = unique([stepTasks{:}]);
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
            for w_ = 1:length(cacheWiring)
                % The caller holds its thread for the whole fetch
                if cacheWiring(w_).readStep == stp && ~isempty(cacheWiring(w_).fetch) && ...
                        ~any(cellfun(@(x) x == cacheWiring(w_).fetch, regionNodes))
                    regionNodes{end+1} = cacheWiring(w_).fetch; %#ok<AGROW>
                end
            end
        end
        % A pooled task keeps one row whose bound covers all its replicas
        tsel_ = ekeyIdx(fcrList(tsel));
        bvec(tsel) = lsn.mult(tsel_) * poolFactor(tsel_);
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
            node = thinkNode{ekey(stepRefTask(i, 1), stepRefTask(i, 2))};
        else
            node = hostStation{stepHost(i, 1), stepHost(i, 2)+1};
        end
    end

    function id = addStep(aidx, hidx, svc, name, blocks, isthink, refKey)
        stepAidx(end+1) = ifempty(aidx, 0); %#ok<AGROW>
        stepHost(end+1, :) = ifempty(hidx, [0, 0]); %#ok<AGROW>
        stepSvc{end+1} = svc; %#ok<AGROW>
        stepName{end+1} = name; %#ok<AGROW>
        stepBlocks(end+1) = blocks; %#ok<AGROW>
        stepIsThink(end+1) = isthink; %#ok<AGROW>
        stepRefTask(end+1, :) = refKey; %#ok<AGROW>
        stepNode{end+1} = []; %#ok<AGROW>
        id = length(stepAidx);
        stepClassOwner(end+1) = id; %#ok<AGROW>
        % Every thread-pool task replica on the stack holds a thread here (a synchronous caller releases it only on reply)
        if isempty(threadStack)
            stepTasks{end+1} = []; %#ok<AGROW>
        else
            stepTasks{end+1} = unique(threadStack(fcrTask(ekeyIdx(threadStack)))); %#ok<AGROW>
        end
    end

    function id = addAuxStep(nodeObj, ownerStep, name, refKey)
        % A step on a Fork, Join or Router node: no station, no service, and
        % no class of its own.
        id = addStep([], [], [], name, false, false, refKey);
        stepNode{id} = nodeObj;
        stepClassOwner(id) = stepClassOwner(ownerStep);
    end

    function [firstStep, replyExits, terminals] = expandEntry(eidx, refKey, trep)
        % Expands the activity subgraph bound to an entry of replica trep of
        % its task, in the call context given by the current entry stack.
        firstStep = [];
        replyExits = zeros(0, 3);
        terminals = zeros(0, 3);

        if ~isempty(entryStack) && any(entryStack(:, 1) == eidx & entryStack(:, 2) == trep)
            line_warning(mfilename, sprintf('Recursive call cycle at entry %s truncated.', lsn.names{eidx}));
            return;
        end
        entryStack(end+1, :) = [eidx, trep];
        threadStack(end+1) = ekey(lsn.parent(eidx), trep);
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

        [firstStep, replyExits, terminals] = expandActivities(boundActs(1), eidx, refKey, trep);

        % Forwarding: with prob p the entry hands off to another entry, which replies directly to the original caller
        fwd = forwardingOf(eidx);
        if ~isempty(fwd) && ~isempty(replyExits)
            ownPorts = replyExits;
            fwdExits = zeros(0, 3);
            pforw = 0.0;
            % Forwarder's thread released at handoff; the forwarded chain expands without it on the thread stack
            fwdThread = threadStack(end);
            threadStack(end) = [];
            fwdTidx = lsn.parent(eidx);
            for f = 1:size(fwd, 1)
                % A forwarding call spreads over the reached callee replicas exactly as a synchronous one does
                reps = targetReplicas(fwdTidx, trep, lsn.parent(fwd(f, 1)));
                p = fwd(f, 2);
                reached = 0;
                for fmrep_ = reps
                    [fFirst, fReplies, fTerms] = expandEntry(fwd(f, 1), refKey, fmrep_);
                    if isempty(fFirst)
                        continue;
                    end
                    reached = reached + 1;
                    for r = 1:size(ownPorts, 1)
                        addRoute(ownPorts(r, 1:2), fFirst, ownPorts(r, 3) * p / length(reps));
                    end
                    fwdExits = [fwdExits; fReplies; fTerms]; %#ok<AGROW>
                end
                pforw = pforw + p * reached / length(reps);
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
        entryStack(end, :) = [];
        threadStack(end) = [];
    end

    function [firstStep, replyExits, terminals] = expandActivities(a0, eidx, refKey, trep)
        % Walks the intra-task activity graph from a0, creating steps and
        % expanding every synchronous call site. Reply exits and terminals
        % are returned as ports, [step, isSignal] rows.
        firstStep = [];
        replyExits = zeros(0, 3);
        terminals = zeros(0, 3);

        tidx = lsn.parent(a0);
        localActs = lsn.actsof{eidx};
        visited = configureDictionary('double','cell'); % aidx -> [entryStep, exitStep, exitIsSignal]
        joinOf = configureDictionary('double','cell');  % join activity -> [joinStep, entryStep]
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
                se = visited{aidx};
                entryStep = se(1);
                exitPort = se(2:3);
                return;
            end
            % The activity bound to an ItemEntry of a CacheTask is the read
            % step: it sits on the Cache node rather than on the processor.
            cacheNode = [];
            if tidx <= length(lsn.iscache) && lsn.iscache(tidx) && full(lsn.graph(eidx, aidx)) > 0
                cacheNode = getCacheNode(tidx, trep);
            end
            [entryStep, exitPort] = makeActivitySteps(aidx, tidx, refKey, trep, cacheNode);
            visited{aidx} = [entryStep, exitPort];

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
                    trigH = addStep(aidx, hostKey(tidx, trep), [], ...
                        sname([lsn.names{aidx}, '_ph2h']), false, false, refKey);
                    trigM = addStep(aidx, hostKey(tidx, trep), [], ...
                        sname([lsn.names{aidx}, '_ph2m']), false, false, refKey);
                    addCacheRoute(entryStep, trigH);
                    addCacheRoute(entryStep, trigM);
                    if hasRetrievalOf(tidx)
                        line_warning(mfilename, sprintf(['Delayed-hit retrieval of %s is not ' ...
                            'represented on a cache read with phase-2 successors.'], ...
                            lsn.names{tidx}));
                    end
                    cacheWiring(end+1) = struct('node', cacheNode, ...
                        'readStep', entryStep, 'hitStep', trigH, 'missStep', trigM, ...
                        'itemproc', lsn.itemproc{eidx}, 'nitems', lsn.nitems(eidx), ...
                        'fetch', [], 'fetchSvc', []); %#ok<AGROW>
                    replyExits(end+1, :) = [trigH, 0, 1.0]; %#ok<AGROW>
                    replyExits(end+1, :) = [trigM, 0, 1.0]; %#ok<AGROW>
                    savedStack = threadStack;
                    threadStack = ekey(tidx, trep);
                    nT0 = size(terminals, 1);
                    for hm = 1:2
                        [sEntry, ~] = walk(succ(hm));
                        if ~isempty(stepNode{sEntry})
                            hmHead = addStep(aidx, hostKey(tidx, trep), [], ...
                                sname(sprintf('%s_ph2b%d', lsn.names{aidx}, hm)), false, false, refKey);
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
                        ph2Exits(end+1, :) = [ph2New(r2, 1:3), refKey]; %#ok<AGROW>
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
                    trig = addStep(aidx, hostKey(tidx, trep), [], ...
                        sname([lsn.names{aidx}, '_ph2t']), false, false, refKey);
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
                    threadStack = ekey(tidx, trep);
                    nT0 = size(terminals, 1);
                    posSucc = succ(full(lsn.graph(aidx, succ)) > 0);
                    target = [];
                    if isAndFork(succ)
                        % Phase 2 opens with an AND-fork: spawn into an immediate head step and fork from there
                        head = addStep(aidx, hostKey(tidx, trep), [], ...
                            sname([lsn.names{aidx}, '_ph2']), false, false, refKey);
                        wireAndFork([head, 0], succ, aidx);
                        target = head;
                    elseif isAndJoinPre(aidx)
                        % Phase 2 at an AND-join branch tail: spawn into an immediate head standing in for this branch at the Join
                        head = addStep(aidx, hostKey(tidx, trep), [], ...
                            sname([lsn.names{aidx}, '_ph2']), false, false, refKey);
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
                        head = addStep(aidx, hostKey(tidx, trep), [], ...
                            sname([lsn.names{aidx}, '_ph2']), false, false, refKey);
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
                        ph2Exits(end+1, :) = [ph2New(r2, 1:3), refKey]; %#ok<AGROW>
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
                    fetch_ = []; fetchSvc_ = [];
                    if hasRetrievalOf(tidx)
                        % The fetch is what the miss branch does, so its demand moves to
                        % the fetch station where concurrent misses coalesce
                        fetch_ = getFetchNode(tidx, trep);
                        fetchSvc_ = stepSvc{mEntry};
                        stepSvc{mEntry} = [];
                        if succ(2) <= length(lsn.callsof) && ~isempty(lsn.callsof{succ(2)})
                            line_warning(mfilename, sprintf(['Calls of miss activity %s stay ' ...
                                'outside the retrieval system, so they are not coalesced ' ...
                                'across concurrent misses.'], lsn.names{succ(2)}));
                        end
                    end
                    cacheWiring(end+1) = struct('node', cacheNode, ...
                        'readStep', entryStep, 'hitStep', hEntry, 'missStep', mEntry, ...
                        'itemproc', lsn.itemproc{eidx}, 'nitems', lsn.nitems(eidx), ...
                        'fetch', fetch_, 'fetchSvc', fetchSvc_); %#ok<AGROW>
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
            forkName_ = sname(['Fork_', lsn.names{aidx}]);
            forkNode = Fork(model, forkName_);
            forkStep = addAuxStep(forkNode, fromPort(1), forkName_, refKey);
            addRoute(fromPort, forkStep, 1.0);
            forkOwnerStack(end+1) = forkStep;
            % Walk a replying branch first so the Join/post-join subgraph is created in its phase-2 context
            repliesFirst_ = arrayfun(@branchReplies, fsucc);
            fsucc = [fsucc(repliesFirst_), fsucc(~repliesFirst_)];
            for b = 1:length(fsucc)
                routerName_ = sname(sprintf('Fork_%s_%d', lsn.names{aidx}, b));
                routerNode = Router(model, routerName_);
                routerStep = addAuxStep(routerNode, forkStep, routerName_, refKey);
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
                js = joinOf{joinAidx};
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
            joinName_ = sname(['Join_', lsn.names{joinAidx}]);
            joinNode = Join(model, joinName_, stepNode{forkOwnerStack(end)});
            joinStep = addAuxStep(joinNode, forkOwnerStack(end), joinName_, refKey);
            addRoute(fromPort, joinStep, 1.0);
            [sEntry, ~] = walk(joinAidx);
            addRoute([joinStep, 0], sEntry, 1.0);
            joinOf{joinAidx} = [joinStep, sEntry];
            joinQuorum(end+1, :) = [joinStep, joinAidx]; %#ok<AGROW>
        end

        function nm = sname(name)
            % A class name carries the replica of the task that owns the step
            nm = uniqueName(suffixed(name, trep));
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

    function [entryStep, exitPort] = makeActivitySteps(aidx, tidx, refKey, trep, cacheNode)
        % One step for the host demand, plus one step per unrolled
        % synchronous call stage.
        hidx = hostKey(tidx, trep);
        if nargin >= 5 && ~isempty(cacheNode)
            % A read step holds no demand and issues no call: the lookup is
            % instantaneous and the work is done on the hit or miss branch.
            entryStep = addStep(aidx, hidx, [], uniqueName(suffixed(lsn.names{aidx}, trep)), false, false, refKey);
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
        entryStep = addStep(aidx, hidx, svc, uniqueName(suffixed(lsn.names{aidx}, trep)), false, false, refKey);
        % cur is the port through which the activity is currently left. A
        % blocking call site is left through its reply signal.
        cur = [entryStep, 0];

        % Activity think time: a delay in series with the host demand, on a shared INF station so the processor is released
        actThinkDist_ = actThinkOf(aidx);
        if ~isempty(actThinkDist_)
            actThinkStep_ = addStep(aidx, hidx, actThinkDist_, ...
                uniqueName(suffixed([lsn.names{aidx}, '_think'], trep)), false, false, refKey);
            stepNode{actThinkStep_} = actThinkStation();
            addRoute(cur, actThinkStep_, 1.0);
            cur = [actThinkStep_, 0];
        end

        % A call blocks the caller's server only when host servers are finite, the task is not a thread pool,
        % multiplicity is finite, and the chain is not open -- see _kb/04-networkstruct.md
        hostBlocks = ~hostIsDelay(hidx(1), hidx(2)+1) && ~fcrTask(tidx) && refKey(1) ~= 0 && ...
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
            % The stage is routed to the callee replicas this caller replica reaches, which share the call mean uniformly
            reps = targetReplicas(tidx, trep, lsn.parent(stage.targetEidx));
            calleeFirsts = [];
            calleeReplies = zeros(0, 3);
            for mrep_ = reps
                [cFirst, cReplies, cTerms] = expandEntry(stage.targetEidx, refKey, mrep_);
                if isempty(cFirst)
                    continue;
                end
                calleeFirsts(end+1) = cFirst; %#ok<AGROW>
                % A callee path that neither replies nor continues still holds a
                % token, so it returns to the caller like a reply would.
                calleeReplies = [calleeReplies; cReplies; cTerms]; %#ok<AGROW>
            end
            if stage.isasync
                threadStack(end+1) = asyncThread;
            end
            if isempty(calleeFirsts)
                continue;   % callee not expandable: drop the call, never block
            end
            share = 1.0 / length(calleeFirsts);

            % Merge step needed when the call may be skipped, another stage follows, or an AND-join branch tail must reach
            % the Join in an ordinary class (a REPLY signal carries no forked-task identity) -- see _kb/04-networkstruct.md
            needsMerge = (stage.prob < 1.0) || (k < length(callStages)) || ...
                ~blocks || isAndJoinPre(aidx);
            if needsMerge
                retName_ = uniqueName(suffixed(sprintf('%s_c%d_ret', lsn.names{aidx}, k), trep));
                nxt = addStep(aidx, hidx, [], retName_, false, false, refKey);
                if ~blocks
                    % Non-blocking return carries no signal, so the merge point can sit on a Router (no server to queue behind)
                    stepNode{nxt} = Router(model, retName_);
                end
            end

            if blocks
                % First mandatory call bound to the service class itself; an intervening class switch would release the server
                if k == 1 && stage.prob >= 1.0 && isequal(cur, [entryStep, 0])
                    blk = entryStep;
                    stepBlocks(blk) = true;
                else
                    blk = addStep(aidx, hidx, [], ...
                        uniqueName(suffixed(sprintf('%s_c%d', lsn.names{aidx}, k), trep)), true, false, refKey);
                    addRoute(cur, blk, stage.prob);
                    if stage.prob < 1.0
                        addRoute(cur, nxt, 1.0 - stage.prob);
                    end
                end
                % Every reached replica replies into the same signal, so the call site blocks once however many replicas it has
                for cf_ = calleeFirsts
                    addRoute([blk, 0], cf_, share);
                end
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
                for cf_ = calleeFirsts
                    addRoute(cur, cf_, stage.prob * share);
                end
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

    function cnode = getCacheNode(tidx, trep)
        % One Cache node per CacheTask replica. The name is suffixed so that it
        % does not collide with the station of the task's processor.
        if isKey(cacheNodeOf, ekey(tidx, trep))
            cnode = cacheNodeOf{ekey(tidx, trep)};
            return;
        end
        cnode = Cache(model, suffixed([lsn.names{tidx}, '_Cache'], trep), lsn.nitems(tidx), ...
            lsn.itemcap{tidx}, lsn.replacestrat(tidx));
        cacheNodeOf{ekey(tidx, trep)} = cnode;
    end

    function [setup_, delayoff_, status_] = functionTimesOf(tidx)
        % Setup and delay-off of a SetupTask. Both are needed: without a
        % delay-off the server never shuts down, so it pays the setup once at
        % most and the pair carries no information
        setup_ = []; delayoff_ = []; status_ = 'none';
        if ~isfield(lsn, 'hassetup') || isempty(lsn.hassetup) || ...
                tidx > numel(lsn.hassetup) || full(lsn.hassetup(tidx)) == 0
            return;
        end
        if ~isfield(lsn, 'setuptime') || tidx > numel(lsn.setuptime)
            return;
        end
        su_ = lsn.setuptime{tidx};
        if ~isa(su_, 'Distribution') || isa(su_, 'Immediate') || su_.getMean() <= GlobalConstants.FineTol
            return;
        end
        do_ = [];
        if isfield(lsn, 'delayofftime') && tidx <= numel(lsn.delayofftime)
            do_ = lsn.delayofftime{tidx};
        end
        if ~isa(do_, 'Distribution')
            status_ = 'nodelayoff';
            return;
        end
        setup_ = su_; delayoff_ = do_; status_ = 'ok';
    end

    function tf = hasRetrievalOf(tidx)
        % True when the CacheTask coalesces concurrent misses of the same item
        tf = isfield(lsn, 'hasretrieval') && ~isempty(lsn.hasretrieval) && ...
            tidx <= length(lsn.hasretrieval) && full(lsn.hasretrieval(tidx)) ~= 0;
    end

    function fnode = getFetchNode(tidx, trep)
        % The retrieval system of a CacheTask is an ordinary queueing network: one
        % PS fetch station per cache replica, as in the LN cache sublayer
        if isKey(fetchNodeOf, ekey(tidx, trep))
            fnode = fetchNodeOf{ekey(tidx, trep)};
            return;
        end
        fnode = Queue(model, suffixed([lsn.names{tidx}, '_Cache_Fetch'], trep), SchedStrategy.PS);
        fetchNodeOf{ekey(tidx, trep)} = fnode;
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

    function applyJoinQuorum(joinNode, joinClass, joinAidx)
        % A join whose quorum equals its branch count already waits for all
        % branches, the default JoinStrategy.STD, so only a genuine quorum
        % k < n switches the node to JoinStrategy.PARTIAL.
        if isempty(joinNode) || isempty(joinClass) || ~isfield(lsn, 'actquorum') || ...
                isempty(lsn.actquorum) || joinAidx > length(lsn.actquorum)
            return;
        end
        quorum = full(lsn.actquorum(joinAidx));
        nbranches = countAndJoinBranches(joinAidx);
        if quorum < 1 || nbranches < 1 || quorum >= nbranches
            return;
        end
        joinNode.setStrategy(joinClass, JoinStrategy.PARTIAL);
        joinNode.setRequired(joinClass, quorum);
    end

    function nb_ = countAndJoinBranches(joinAidx)
        % Branch tails feeding an AND-join, i.e. its PRE_AND predecessors.
        nb_ = 0;
        for p_ = find(lsn.graph(:, joinAidx) ~= 0)'
            if p_ ~= joinAidx && isAndJoinPre(p_)
                nb_ = nb_ + 1;
            end
        end
    end

    function ad_ = actThinkOf(aidx)
        % Think time of an activity, empty when it has none. It is a delay in
        % series with the host demand, held at the activity's own task (the
        % thread is kept) but with the host processor released, as in lqns.
        ad_ = [];
        if ~isfield(lsn, 'actthink') || aidx > length(lsn.actthink)
            return;
        end
        at_ = lsn.actthink{aidx};
        if isa(at_, 'Distribution') && ~isa(at_, 'Immediate') && at_.getMean() > GlobalConstants.FineTol
            ad_ = at_;
        end
    end

    function nd_ = actThinkStation()
        % Single INF station shared by every activity think time.
        if isempty(actThinkNode)
            actThinkNode = Delay(model, 'ActivityThink');
        end
        nd_ = actThinkNode;
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

    function n_ = nrep(idx)
        % Replicas materialised for a host or task index
        if materialize
            n_ = replRaw(idx);
        else
            n_ = 1;
        end
    end

    function f_ = poolFactor(idx)
        % Capacity multiplier carried by a pooled element, 1 when materialised
        if materialize
            f_ = 1;
        else
            f_ = replRaw(idx);
        end
    end

    function f_ = fanOutOf(aTidx, bTidx, rb)
        % Callee replicas reached by one caller replica. An unset fan-out is the
        % smallest value consistent with repl(a)*fanout = repl(b)*fanin.
        if rb <= 1
            f_ = 1;
            return;
        end
        f_ = 0;
        if isfield(lsn, 'fanout') && ~isempty(lsn.fanout) && ...
                aTidx <= size(lsn.fanout, 1) && bTidx <= size(lsn.fanout, 2)
            f_ = full(lsn.fanout(aTidx, bTidx));
        end
        if f_ <= 0
            ra = replRaw(aTidx);
            if rb > ra
                f_ = max(1, floor(rb / ra));
            else
                f_ = 1;
            end
        end
        f_ = min(max(1, round(f_)), rb);
    end

    function reps = targetReplicas(aTidx, aRep, bTidx)
        % Replicas of the callee reached by replica aRep of the caller
        rb = nrep(bTidx);
        if rb <= 1
            reps = 0;
            return;
        end
        f_ = fanOutOf(aTidx, bTidx, rb);
        reps = mod(aRep * f_ + (0:f_-1), rb);
    end

    function k_ = hostKey(tidx_, trep_)
        % Station of the processor replica running replica trep of a task
        hidx_ = lsn.parent(tidx_);
        k_ = [hidx_, mod(trep_, nrep(hidx_))];
    end

    function k_ = ekey(idx, rep)
        % Composite (element, 0-based replica) key for maps and stacks
        k_ = (idx - 1) * RKEY + rep + 1;
    end

    function idx = ekeyIdx(k_)
        % Element index of a composite key
        idx = floor((k_ - 1) / RKEY) + 1;
    end

    function nm = uniqueName(base)
        % A subgraph copied per call site or per replica repeats its names, which
        % Network rejects, so a repeat is disambiguated by occurrence
        if isKey(usedNames, base)
            n_ = usedNames(base) + 1;
            usedNames(base) = n_;
            nm = sprintf('%s_d%d', base, n_);
        else
            usedNames(base) = 1;
            nm = base;
        end
    end

    function nm = suffixed(name, rep)
        % Replica 1 keeps the plain name. The suffix stays inside [A-Za-z0-9_]
        % because a class name is a JSON object key, and jsondecode mangles any
        % key that is not a valid identifier -- see _kb/11-conventions-and-gotchas.md
        if rep == 0
            nm = name;
        else
            nm = sprintf('%s_r%d', name, rep + 1);
        end
    end

    function n_ = replInstantiations()
        % Copies of the task step graphs a materialised expansion would create
        rimemo_ = configureDictionary('double', 'double');
        n_ = 0;
        riseeds_ = unique([refTaskIndices(:).', arrayfun(@(e) lsn.parent(e), openEntries)]);
        for rit_ = riseeds_
            n_ = n_ + replRaw(rit_) * taskCost(rit_, []);
        end

        function ric_ = taskCost(ritidx_, ristack_)
            if any(ristack_ == ritidx_)
                ric_ = 1;   % recursive cycle: truncated anyway
                return;
            end
            if isKey(rimemo_, ritidx_)
                ric_ = rimemo_(ritidx_);
                return;
            end
            rideeper_ = [ristack_, ritidx_];
            ric_ = 1;
            for rieidx_ = lsn.entriesof{ritidx_}
                ritargets_ = [];
                for ria_ = lsn.actsof{rieidx_}
                    if ria_ > length(lsn.callsof) || isempty(lsn.callsof{ria_})
                        continue;
                    end
                    for ricidx_ = lsn.callsof{ria_}
                        rict_ = full(lsn.calltype(ricidx_));
                        if rict_ == CallType.SYNC || rict_ == CallType.ASYNC
                            ritargets_(end+1) = lsn.callpair(ricidx_, 2); %#ok<AGROW>
                        end
                    end
                end
                % A forwarded entry is expanded per replica just as a call is
                rifwd_ = forwardingOf(rieidx_);
                if ~isempty(rifwd_)
                    ritargets_ = [ritargets_, rifwd_(:, 1).']; %#ok<AGROW>
                end
                for rite_ = ritargets_
                    rib_ = lsn.parent(rite_);
                    ric_ = ric_ + fanOutOf(ritidx_, rib_, replRaw(rib_)) * taskCost(rib_, rideeper_);
                end
            end
            rimemo_(ritidx_) = ric_;
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
