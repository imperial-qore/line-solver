function  buildLayersRecursive(self, idxSet, callers, ishostlayer, flat)
% IDXSET is the server element of this layer, a scalar under 'srvn' layering
% and the whole host+task set under 'flat' layering; see _kb/06-solver-catalog.md
if nargin < 5
    flat = false;
end
lqn = self.lqn;
idx = idxSet(1); % layer key: model name, ensemble slot and update-map column
jobPosKey = zeros(lqn.nidx,1);
curClassKey = cell(lqn.nidx,1);
curStationKey = cell(lqn.nidx,1);
% Fan-out check: see _kb/06-solver-catalog.md (LN section) for rationale
rawReplicas = lqn.repl(idx);
reduceFanout = false;
if ~flat && rawReplicas > 1 && ~isempty(callers)
    if ~ishostlayer && isfield(lqn, 'fanout') && ~isempty(lqn.fanout)
        reduceFanout = true;
        for c = callers(:)'
            if lqn.fanout(c, idx) < rawReplicas
                reduceFanout = false;
                break;
            end
        end
    elseif ishostlayer
        reduceFanout = true;
        for c = callers(:)'
            if lqn.repl(c) ~= rawReplicas
                reduceFanout = false;
                break;
            end
        end
    end
end
if reduceFanout
    nreplicas = 1;
    if ~ishostlayer
        self.singleReplicaTasks(end+1) = idx;
    end
elseif flat
    nreplicas = 1; % flat layering rejects replicated elements upstream
else
    nreplicas = rawReplicas;
end
%mult = lqn.mult;
mult = lqn.maxmult; % this removes spare capacity that cannot be used
lqn.mult = mult;
callservtproc = self.callservtproc;
if flat
    model = Network([lqn.hashnames{idx},'.Flat']);
else
    model = Network(lqn.hashnames{idx});
end
model.setChecks(false); % fast mode
model.attribute = struct('hosts',[],'tasks',[],'entries',[],'activities',[],'calls',[],'serverIdx',0);
model.attribute.serverIdxOf = NaN(lqn.nidx,1); % LQN server element -> its station index
if flat | ishostlayer | any(any(lqn.issynccaller(callers, entriesOfSet()))) | any(any(lqn.isasynccaller(callers, entriesOfSet()))) %#ok<OR2>
    clientDelay = Delay(model, 'Clients');
    model.attribute.clientIdx = 1;
    model.attribute.serverIdx = 2;
    model.attribute.sourceIdx = NaN;
else
    model.attribute.serverIdx = 1;
    model.attribute.clientIdx = NaN;
    model.attribute.sourceIdx = NaN;
end
% One station (times its replicas) per server element of the layer
srvStation = cell(lqn.nidx,1);
% (a setup no longer changes how a layer is built: see lqn_setup_charge)
for sidx = idxSet(:)'
    sishost = sidx <= lqn.nhosts;
    srvStation{sidx} = cell(1,nreplicas);
    for m=1:nreplicas
        if m == 1
            srvStation{sidx}{m} = Queue(model,lqn.hashnames{sidx}, lqn.sched(sidx));
        else
            srvStation{sidx}{m} = Queue(model,[lqn.hashnames{sidx},'.',num2str(m)], lqn.sched(sidx));
        end
        srvStation{sidx}{m}.setNumberOfServers(mult(sidx));
        srvStation{sidx}{m}.attribute.ishost = sishost;
        srvStation{sidx}{m}.attribute.idx = sidx;
        % see _kb/06-solver-catalog.md (LN section) for rationale
        srvStation{sidx}{m}.setImmediateFeedback(true);
    end
    model.attribute.serverIdxOf(sidx) = model.getNodeIndex(srvStation{sidx}{1});
end
serverStation = srvStation{idx}; % the layer's own server, sole server under 'srvn'

iscachelayer = ~flat && all(lqn.iscache(callers)) && ishostlayer;
hasRetrievalCache = false;
retrievalStation = [];
if iscachelayer
    cacheNode = Cache(model, lqn.hashnames{callers}, lqn.nitems(callers), lqn.itemcap{callers}, lqn.replacestrat(callers));
    % Delayed-hit retrieval: a dedicated fetch station in the cache sublayer so the
    % closed AMVA (da_cacheqn_retrieval) captures the finite-population coalescing.
    hasRetrievalCache = isfield(lqn,'hasretrieval') && any(lqn.hasretrieval(callers));
    if hasRetrievalCache
        retrievalStation = Queue(model, [lqn.hashnames{callers},'.Fetch'], SchedStrategy.PS);
    end
end
retrievalWiring = [];

actsInCaller = [lqn.actsof{callers}];
isPostAndAct = full(lqn.actposttype)==ActivityPrecedenceType.POST_AND;
isPreAndAct = full(lqn.actpretype)==ActivityPrecedenceType.PRE_AND;
hasfork = any(intersect(find(isPostAndAct),actsInCaller));

maxfanout = 1; % maximum output parallelism level of fork nodes
for aidx = actsInCaller(:)'
    successors = find(lqn.graph(aidx,:));
    if any(isPostAndAct(successors))
        maxfanout = max(maxfanout, sum(isPostAndAct(successors)));
    end
end

if hasfork
    forkNode = Fork(model, 'Fork_PostAnd');
    for f=1:maxfanout
        forkOutputRouter{f} = Router(model, ['Fork_PostAnd_',num2str(f)]);
    end
    forkClassStack = []; % stack with the entry class at the visited forks, the last visited is end of the list.
end

isPreAndAct = full(lqn.actpretype)==ActivityPrecedenceType.PRE_AND;
hasjoin = any(isPreAndAct(actsInCaller));
if hasjoin
    joinNode = Join(model, 'Join_PreAnd', forkNode);
end

aidxClass = cell(1, lqn.nidx);
aidxThinkClass = cell(1, lqn.nidx); % auxiliary classes for activity think-time
cidxClass = cell(1,0);
cidxAuxClass = cell(1,0);

% Routed call groups: a group is ONE dispatch with n destinations, so its
% members share a call class and a dispatch class that holds the job at the
% client while the target is picked. See _kb/06-solver-catalog.md (LN section).
[cgroupOfCidx, cgroupMembers, cgroupStrategy] = callGroupsByCidx(lqn);
grpDispatchClass = cell(1, numel(cgroupMembers));
grpCallClass = cell(1, numel(cgroupMembers));
grpRouter = cell(1, numel(cgroupMembers));
groupRouted = false(1, numel(cgroupMembers));
routedGroupSites = zeros(0,3); % [router node index, dispatch class index, group index]

self.servt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % server classes to update
self.thinkt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % client classes to update
self.actthinkt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % activity think-time classes to update
self.arvproc_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % classes to update in the next iteration for asynch calls
self.call_classes_updmap{idx} = zeros(0,4); % [modelidx, callidx, node, class] % calls classes to update in the next iteration (includes calls in client classes)
self.route_prob_updmap{idx} = zeros(0,7); % [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto] % routing probabilities to update in the next iteration

% Station indices of the layer's servers, kept apart from attribute.hosts /
% attribute.tasks whose rows are [class index, LQN element] pairs
model.attribute.hostStations = [];
model.attribute.taskStations = [];
for sidx = idxSet(:)'
    if sidx <= lqn.nhosts
        model.attribute.hostStations(end+1) = model.attribute.serverIdxOf(sidx);
        if ~flat
            model.attribute.hosts(end+1,:) = [NaN, model.attribute.serverIdxOf(sidx) ];
        end
    else
        model.attribute.taskStations(end+1) = model.attribute.serverIdxOf(sidx);
        if ~flat
            model.attribute.tasks(end+1,:) = [NaN, model.attribute.serverIdxOf(sidx) ];
        end
    end
end

hasSource = false; % flag whether a source is needed
openClasses = [];
entryOpenClasses = []; % track entry-level open arrivals
% first pass: create the classes
for tidx_caller = callers
    % For host layers, check if the task has any entries with sync/async callers
    % or has open arrivals, OR if any entry is a forwarding target.
    hasDirectCallers = false;
    isForwardingTarget = false;
    hasOpenArrival = false;
    if hostIsServer(tidx_caller)
        % Check if the task is a reference task (always create closed class)
        if lqn.isref(tidx_caller)
            hasDirectCallers = true;
        else
            % Check if any entry of this task has sync or async callers
            for eidx = lqn.entriesof{tidx_caller}
                if any(full(lqn.issynccaller(:, eidx))) || any(full(lqn.isasynccaller(:, eidx)))
                    hasDirectCallers = true;
                end
                % Also check for open arrivals on this entry
                if isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && ...
                   iscell(lqn.arrival) && eidx <= length(lqn.arrival) && ...
                   ~isempty(lqn.arrival{eidx})
                    hasOpenArrival = true;
                end
                % Check if this entry is a forwarding target
                for cidx = 1:lqn.ncalls
                    if full(lqn.calltype(cidx)) == CallType.FWD && full(lqn.callpair(cidx, 2)) == eidx
                        isForwardingTarget = true;
                        break;
                    end
                end
            end
        end
    end
    isLayerClient = any(any(lqn.issynccaller(tidx_caller, entriesOfSet())));
    % A task reached ONLY by an entry arrival has no task layer, because no task
    % calls it, so updateThinkTimes never gives its caller class a surrogate delay
    % and the class cycles against an Immediate one. Adding an open stream on top
    % of that unthrottled chain is what saturated lqn_open_arrival: 0.68 of the
    % processor under 'srvn' and 1.000 under 'srvn.ph', against 0.32 from lqns,
    % lqsim and LDES alike. The chain is the better of the two representations
    % here -- it is the one that honours the thread pool, which is where the
    % requests actually queue -- so the stream is dropped and the chain is closed
    % on the known arrival rate instead, exactly as a forwarding target is.
    openArrivalOnly = hasOpenArrival && ~hasDirectCallers && ~isForwardingTarget ...
        && ~isLayerClient && ~lqn.isref(tidx_caller);
    if (hostIsServer(tidx_caller) && (hasDirectCallers || hasOpenArrival || isForwardingTarget)) | isLayerClient %#ok<OR2> % if it is only an asynch caller the closed classes are not needed
        if self.njobs(tidx_caller,idx) == 0
            % for each entry of the calling task
            % determine job population
            % this block matches the corresponding calculations in
            % updateThinkTimes
            % Use single-replica njobs if this layer or the caller is in single-replica mode
            callerIsSingleReplica = reduceFanout || any(self.singleReplicaTasks == tidx_caller);
            if callerIsSingleReplica
                njobs = mult(tidx_caller);
            else
                njobs = mult(tidx_caller)*lqn.repl(tidx_caller);
            end
            if isinf(njobs)
                callers_of_tidx_caller = find(lqn.taskgraph(:,tidx_caller));
                njobs = sum(mult(callers_of_tidx_caller)); %#ok<FNDSB>
                if isinf(njobs)
                    % if also the callers of tidx_caller are inf servers, then use
                    % an heuristic
                    njobs = min(sum(mult(isfinite(mult)) .* lqn.repl(isfinite(mult))),1000); % Python parity: cap at 1000
                end
            end
            self.njobs(tidx_caller,idx) = njobs;
        else
            njobs = self.njobs(tidx_caller,idx);
        end
        caller_name = lqn.hashnames{tidx_caller};
        aidxClass{tidx_caller} = ClosedClass(model, caller_name, njobs, clientDelay);
        clientDelay.setService(aidxClass{tidx_caller}, Disabled.getInstance());
        setAllServers(aidxClass{tidx_caller}, Disabled.getInstance());
        aidxClass{tidx_caller}.completes = false;
        aidxClass{tidx_caller}.setReferenceClass(true); % renormalize residence times using the visits to the task
        aidxClass{tidx_caller}.attribute = [LayeredNetworkElement.TASK, tidx_caller];
        model.attribute.tasks(end+1,:) = [aidxClass{tidx_caller}.index, tidx_caller];
        if lqn.isref(tidx_caller)
            clientDelay.setService(aidxClass{tidx_caller}, self.thinkproc{tidx_caller});
        else
            % a served task's declared think time is not a per-request delay,
            % so the seed carries none either; updateThinkTimes replaces this
            % with the surrogate delay from the first iteration on
            clientDelay.setService(aidxClass{tidx_caller}, Immediate.getInstance());
            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, tidx_caller, 1, aidxClass{tidx_caller}.index];
        end
        for eidx = lqn.entriesof{tidx_caller}
            % create a class
            aidxClass{eidx} = ClosedClass(model, lqn.hashnames{eidx}, 0, clientDelay);
            clientDelay.setService(aidxClass{eidx}, Disabled.getInstance());
            setAllServers(aidxClass{eidx}, Disabled.getInstance());
            aidxClass{eidx}.completes = false;
            aidxClass{eidx}.attribute = [LayeredNetworkElement.ENTRY, eidx];
            model.attribute.entries(end+1,:) = [aidxClass{eidx}.index, eidx];
            [singleton, javasingleton] = Immediate.getInstance();
            if isempty(model.obj)
                clientDelay.setService(aidxClass{eidx}, singleton);
            else
                clientDelay.setService(aidxClass{eidx}, javasingleton);
            end

            % Check for open arrival distribution on this entry
            if ~openArrivalOnly && isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && ...
               iscell(lqn.arrival) && eidx <= length(lqn.arrival) && ...
               ~isempty(lqn.arrival{eidx})

                if ~hasSource
                    hasSource = true;
                    model.attribute.sourceIdx = length(model.nodes)+1;
                    sourceStation = Source(model,'Source');
                    sinkStation = Sink(model,'Sink');
                end

                % Create open class for this entry
                openClassForEntry = OpenClass(model, [lqn.hashnames{eidx}, '_Open'], 0);
                sourceStation.setArrival(openClassForEntry, lqn.arrival{eidx});
                clientDelay.setService(openClassForEntry, Disabled.getInstance());

                % Use bound activity's service time (entries themselves have Immediate service)
                % Find activities bound to this entry via graph
                bound_act_indices = find(lqn.graph(eidx,:) > 0);
                if ~isempty(bound_act_indices)
                    % Use first bound activity's service time
                    bound_aidx = bound_act_indices(1);
                    setAllServers(openClassForEntry, self.servtproc{bound_aidx});
                else
                    % Fallback to entry service (should not happen in well-formed models)
                    setAllServers(openClassForEntry, self.servtproc{eidx});
                end

                % Track for routing setup later: [class_index, entry_index]
                entryOpenClasses(end+1,:) = [openClassForEntry.index, eidx];

                % Track: Use negative entry index to distinguish from call arrivals
                self.arvproc_classes_updmap{idx}(end+1,:) = [idx, -eidx, ...
                    model.getNodeIndex(sourceStation), openClassForEntry.index];

                openClassForEntry.completes = false;
                openClassForEntry.attribute = [LayeredNetworkElement.ENTRY, eidx];
            end
        end
    end

    % for each activity of the calling task
    for aidx = lqn.actsof{tidx_caller}
        if hostIsServer(tidx_caller) | any(any(lqn.issynccaller(tidx_caller, entriesOfSet()))) %#ok<OR2>
            % create a class
            aidxClass{aidx} = ClosedClass(model, lqn.hashnames{aidx}, 0, clientDelay);
            clientDelay.setService(aidxClass{aidx}, Disabled.getInstance());
            setAllServers(aidxClass{aidx}, Disabled.getInstance());
            aidxClass{aidx}.completes = false;
            aidxClass{aidx}.attribute = [LayeredNetworkElement.ACTIVITY, aidx];
            model.attribute.activities(end+1,:) = [aidxClass{aidx}.index, aidx];
            hidx = lqn.parent(lqn.parent(aidx)); % index of host processor
            if isempty(serversFor(hidx))
                % set the host demand for the activity
                clientDelay.setService(aidxClass{aidx}, self.servtproc{aidx});
            end
            if lqn.sched(tidx_caller)~=SchedStrategy.REF % in 'ref' case the service activity is constant
                % updmap(end+1,:) = [idx, aidx, 1, idxClass{aidx}.index];
            end
            if iscachelayer && full(lqn.graph(eidx,aidx))
                clientDelay.setService(aidxClass{aidx}, self.servtproc{aidx});
            end

            % Auxiliary class carrying the activity think time, host layer only.
            % see _kb/06-solver-catalog.md (LN section) for rationale
            if ~isempty(lqn.actthink{aidx}) && lqn.actthink{aidx}.getMean() > GlobalConstants.FineTol ...
                    && ~isempty(serversFor(hidx))
                aidxThinkClass{aidx} = ClosedClass(model, [lqn.hashnames{aidx},'.Think'], 0, clientDelay);
                aidxThinkClass{aidx}.completes = false;
                aidxThinkClass{aidx}.attribute = [LayeredNetworkElement.ACTIVITY, aidx];
                clientDelay.setService(aidxThinkClass{aidx}, lqn.actthink{aidx});
                setAllServers(aidxThinkClass{aidx}, Disabled.getInstance());
                self.actthinkt_classes_updmap{idx}(end+1,:) = [idx, aidx, 1, aidxThinkClass{aidx}.index];
            end
        end
        % add a class for each outgoing call from this activity
        for cidx = lqn.callsof{aidx}
            callmean(cidx) = lqn.callproc{cidx}.getMean;
            switch lqn.calltype(cidx)
                case CallType.ASYNC
                    if ~isempty(serversFor(lqn.parent(lqn.callpair(cidx,2)))) % add only if the target is a server here
                        if ~hasSource % we need to add source and sink to the model
                            hasSource = true;
                            model.attribute.sourceIdx = length(model.nodes)+1;
                            sourceStation = Source(model,'Source');
                            sinkStation = Sink(model,'Sink');
                        end
                        cidxClass{cidx} = OpenClass(model, lqn.callhashnames{cidx}, 0);
                        sourceStation.setArrival(cidxClass{cidx}, Immediate.getInstance());
                        clientDelay.setService(cidxClass{cidx}, Disabled.getInstance());
                        setAllServers(cidxClass{cidx}, Immediate.getInstance());
                        openClasses(end+1,:) = [cidxClass{cidx}.index, callmean(cidx), cidx];
                        model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                        cidxClass{cidx}.completes = false;
                        cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        seedCallService(cidx);
                    end
                case CallType.SYNC
                    gid = cgroupOfCidx(cidx);
                    if gid > 0
                        % One dispatch with n destinations. The members share the
                        % dispatch class, which is the one the strategy routes and
                        % the one that visits the targets: the hop must NOT switch
                        % class, because a state-dependent routing function is
                        % evaluated at zero off the class diagonal (sub_jsq). The
                        % class switch goes on the return arc instead, into a group
                        % class the job continues in.
                        if isempty(grpDispatchClass{gid})
                            % The strategy is a property of a NODE, and it routes
                            % over that node's links, not over the arcs of one
                            % class. A dedicated router whose only links are the
                            % group's targets is therefore the only place where
                            % the choice is exactly the group's
                            grpRouter{gid} = Router(model, [lqn.hashnames{aidx},'.Dispatch',num2str(gid),'.Router']);
                            grpDispatchClass{gid} = ClosedClass(model, [lqn.hashnames{aidx},'.Dispatch',num2str(gid)], 0, clientDelay);
                            grpDispatchClass{gid}.completes = false;
                            grpDispatchClass{gid}.attribute = [LayeredNetworkElement.CALL, cidx];
                            clientDelay.setService(grpDispatchClass{gid}, Immediate.getInstance());
                            setAllServers(grpDispatchClass{gid}, Disabled.getInstance());
                            grpCallClass{gid} = ClosedClass(model, [lqn.callhashnames{cidx},'.Group',num2str(gid)], 0, clientDelay);
                            clientDelay.setService(grpCallClass{gid}, Immediate.getInstance());
                            setAllServers(grpCallClass{gid}, Disabled.getInstance());
                            grpCallClass{gid}.completes = false;
                            grpCallClass{gid}.attribute = [LayeredNetworkElement.CALL, cidx];
                            routedGroupSites(end+1,:) = [model.getNodeIndex(grpRouter{gid}), grpDispatchClass{gid}.index, gid]; %#ok<AGROW>
                        end
                        cidxClass{cidx} = grpDispatchClass{gid};
                    else
                        cidxClass{cidx} = ClosedClass(model, lqn.callhashnames{cidx}, 0, clientDelay);
                        clientDelay.setService(cidxClass{cidx}, Disabled.getInstance());
                        setAllServers(cidxClass{cidx}, Disabled.getInstance());
                        cidxClass{cidx}.completes = false;
                        cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                    end
                    model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                    seedCallService(cidx);
            end

            % an Aux class is needed whenever the call does not happen exactly
            % once per activity execution, which is a property of callmean alone.
            % A group member's mean is the 1/n share the dispatch already carries,
            % so an Aux there would charge the skip path twice
            if callmean(cidx) ~= 1 && cgroupOfCidx(cidx) == 0
                switch lqn.calltype(cidx)
                    case CallType.SYNC
                        cidxAuxClass{cidx} = ClosedClass(model, [lqn.callhashnames{cidx},'.Aux'], 0, clientDelay);
                        cidxAuxClass{cidx}.completes = false;
                        cidxAuxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        clientDelay.setService(cidxAuxClass{cidx}, Immediate.getInstance());
                        setAllServers(cidxAuxClass{cidx}, Disabled.getInstance());
                end
            end

            % see _kb/06-solver-catalog.md (LN section) for rationale
        end
    end
end

% see _kb/06-solver-catalog.md (LN section) for rationale
if hasSource
    nClasses = model.getNumberOfClasses();
    for k = 1:nClasses
        if k > length(sourceStation.input.sourceClasses) || isempty(sourceStation.input.sourceClasses{k})
            sourceStation.input.sourceClasses{k} = {[], ServiceStrategy.LI, Disabled.getInstance()};
        end
        if k > length(sourceStation.arrivalProcess) || isempty(sourceStation.arrivalProcess{k})
            sourceStation.arrivalProcess{k} = Disabled.getInstance();
        end
    end
end

% The fork-join transform mints its own Source/Sink pair, detaching the open
% stream already routed through this one: see _kb/06-solver-catalog.md (LN section)
if hasSource && hasfork
    line_error(mfilename, ['SolverLN: layer ''%s'' carries both an AND fork and an ' ...
        'open stream (an async call or an entry arrival); the fork-join transform ' ...
        'needs a Source of its own'], model.getName());
end

P = model.initRoutingMatrix;
if hasSource
    for o = 1:size(openClasses,1)
        oidx = openClasses(o,1);
        callmean_o = openClasses(o,2);
        p = 1 / callmean_o; % divide by mean number of calls, they go to a server at random
        cidx = openClasses(o,3); % 3 = source
        tgtStation = serversFor(lqn.parent(lqn.callpair(cidx,2)));
        if callmean_o < 1
            % fewer than one call per arrival: a single Bernoulli pass, the
            % geometric loop below would need a negative repeat probability
            P{model.classes{oidx}, model.classes{oidx}}(sourceStation,sinkStation) = 1 - callmean_o;
            for m=1:length(tgtStation)
                P{model.classes{oidx}, model.classes{oidx}}(sourceStation,tgtStation{m}) = callmean_o/length(tgtStation);
                P{model.classes{oidx}, model.classes{oidx}}(tgtStation{m},sinkStation) = 1;
            end
        else
            for m=1:length(tgtStation)
                P{model.classes{oidx}, model.classes{oidx}}(sourceStation,tgtStation{m}) = 1/length(tgtStation);
                for n=1:length(tgtStation)
                    P{model.classes{oidx}, model.classes{oidx}}(tgtStation{m},tgtStation{n}) = (1-p)/length(tgtStation);
                end
                P{model.classes{oidx}, model.classes{oidx}}(tgtStation{m},sinkStation) = p;
            end
        end
        self.arvproc_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(sourceStation), oidx];
        for m=1:length(tgtStation)
            self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(tgtStation{m}), oidx];
        end
    end
end

%% job positions are encoded as follows: 1=client, 2=any of the nreplicas server stations, 3=cache node, 4=fork node, 5=join node
atClient = 1;
atServer = 2;
atCache = 3;

jobPos = atClient; % start at client
% second pass: setup the routing out of entries
for tidx_caller = callers
    % Use same condition as first pass - only process if closed class was created
    hasDirectCallers = false;
    isForwardingTarget = false;
    if hostIsServer(tidx_caller)
        if lqn.isref(tidx_caller)
            hasDirectCallers = true;
        else
            for eidx_check = lqn.entriesof{tidx_caller}
                if any(full(lqn.issynccaller(:, eidx_check))) || any(full(lqn.isasynccaller(:, eidx_check)))
                    hasDirectCallers = true;
                    break;
                end
                if isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && ...
                   iscell(lqn.arrival) && eidx_check <= length(lqn.arrival) && ...
                   ~isempty(lqn.arrival{eidx_check})
                    hasDirectCallers = true;
                    break;
                end
                % Check if this entry is a forwarding target
                for cidx_fwd = 1:lqn.ncalls
                    if full(lqn.calltype(cidx_fwd)) == CallType.FWD && full(lqn.callpair(cidx_fwd, 2)) == eidx_check
                        isForwardingTarget = true;
                        break;
                    end
                end
            end
        end
    end
    if (hostIsServer(tidx_caller) && (hasDirectCallers || isForwardingTarget)) | any(any(lqn.issynccaller(tidx_caller, entriesOfSet()))) %#ok<OR2>
        % for each entry of the calling task
        ncaller_entries = length(lqn.entriesof{tidx_caller});
        for eidx = lqn.entriesof{tidx_caller}
            aidxClass_eidx = aidxClass{eidx};
            aidxClass_tidx_caller = aidxClass{tidx_caller};
            % initialize the probability to select an entry to be identical
            P{aidxClass_tidx_caller, aidxClass_eidx}(clientDelay, clientDelay) = 1 / ncaller_entries;
            if ncaller_entries > 1
                % at successive iterations make sure to replace this with throughput ratio
                self.route_prob_updmap{idx}(end+1,:) = [idx, tidx_caller, eidx, 1, 1, aidxClass_tidx_caller.index, aidxClass_eidx.index];
            end
            P = recurActGraph(P, tidx_caller, eidx, aidxClass_eidx, jobPos, {});
        end
    end
end

% Setup routing for entry-level open arrivals (AFTER recurActGraph to avoid being overwritten)
if hasSource && ~isempty(entryOpenClasses)
    for e = 1:size(entryOpenClasses,1)
        eoidx = entryOpenClasses(e,1); % class index
        openClass = model.classes{eoidx};

        % Explicitly set routing: ONLY Source → Server → Sink
        % Zero out all routing for this class first
        for node1 = 1:length(model.nodes)
            for node2 = 1:length(model.nodes)
                P{openClass, openClass}(node1, node2) = 0;
            end
        end

        % Now set the correct routing
        eoStation = entryOpenStation(entryOpenClasses(e,2));
        for m=1:length(eoStation)
            % Route: Source -> station of the entry -> Sink
            P{openClass, openClass}(sourceStation,eoStation{m}) = 1/length(eoStation);
            P{openClass, openClass}(eoStation{m},sinkStation) = 1.0;
        end
    end
end

% see _kb/06-solver-catalog.md (LN section) for rationale
if hasSource
    nClasses = model.getNumberOfClasses();
    for k = 1:nClasses
        if k > length(sourceStation.input.sourceClasses) || isempty(sourceStation.input.sourceClasses{k})
            sourceStation.input.sourceClasses{k} = {[], ServiceStrategy.LI, Disabled.getInstance()};
        end
        if k > length(sourceStation.arrivalProcess) || isempty(sourceStation.arrivalProcess{k})
            sourceStation.arrivalProcess{k} = Disabled.getInstance();
        end
    end
end

% Delayed-hit retrieval cache wiring (EXPERIMENTAL).
% see _kb/06-solver-catalog.md (LN section) for rationale
if hasRetrievalCache && ~isempty(retrievalWiring)
    retrievalStation.setService(retrievalWiring.readClass, self.servtproc{retrievalWiring.missaidx});
    P{retrievalWiring.readClass, retrievalWiring.readClass}(cacheNode, retrievalStation) = 1.0;
    P{retrievalWiring.readClass, retrievalWiring.readClass}(retrievalStation, cacheNode) = 1.0;
    cacheNode.setRetrievalSystem(retrievalWiring.readClass, retrievalWiring.missClass, retrievalStation);
    % see _kb/06-solver-catalog.md (LN section) for rationale
    rcvals = unique(cacheNode.server.retrievalClasses(:));
    rcvals = rcvals(rcvals > 0);
    for rci = rcvals(:)'
        model.classes{rci}.attribute = [-1, -1];
        model.classes{rci}.completes = false;
    end
    % see _kb/06-solver-catalog.md (LN section) for rationale
    hc = cacheNode.server.hitClass; mc = cacheNode.server.missClass;
    Lhm = max(numel(hc), numel(mc));
    if numel(hc) < Lhm, hc(numel(hc)+1:Lhm) = 0; cacheNode.server.hitClass = hc; end
    if numel(mc) < Lhm, mc(numel(mc)+1:Lhm) = 0; cacheNode.server.missClass = mc; end
    % see _kb/06-solver-catalog.md (LN section) for rationale
    KfullSvc = model.getNumberOfClasses;
    for istp = 1:numel(model.stations)
        stp = model.stations{istp};
        if isprop(stp,'server') && ~strcmpi(class(stp.server),'ServiceTunnel')
            for rp = 1:KfullSvc
                if numel(stp.server.serviceProcess) < rp || isempty(stp.server.serviceProcess{rp})
                    stp.setService(model.classes{rp}, Disabled.getInstance());
                end
            end
        end
    end
    % see _kb/06-solver-catalog.md (LN section) for rationale
    Kfull = model.getNumberOfClasses;
    Icur = model.getNumberOfNodes;
    if isa(P, 'RoutingMatrix')
        Pc = P.getCell();
    else
        Pc = P;
    end
    if size(Pc, 1) < Kfull
        Pnew = cell(Kfull, Kfull);
        [ro, co] = size(Pc);
        Pnew(1:ro, 1:co) = Pc;
        for rr = 1:Kfull
            for ss = 1:Kfull
                if isempty(Pnew{rr, ss})
                    Pnew{rr, ss} = zeros(Icur, Icur);
                end
            end
        end
        Pc = Pnew;
    end
    P = Pc;
end

if flat
    % link() installs RAND routing for every (node, class) pair left without
    % an outgoing arc. With one station per server in a single layer those
    % spurious uniform arcs let a class wander to stations it never visits,
    % trapping the flow in a sub-cycle and leaving the reference class with
    % zero visits. Disable them, as LQN2QN does for its signal classes.
    if isa(P, 'RoutingMatrix')
        Plink = P.getCell();
    else
        Plink = P;
    end
    for inode = 1:length(model.nodes)
        if isa(model.nodes{inode},'Sink')
            continue
        end
        for rcls = 1:model.getNumberOfClasses
            outflow = 0;
            for scls = 1:model.getNumberOfClasses
                if rcls <= size(Plink,1) && scls <= size(Plink,2) && ~isempty(Plink{rcls,scls})
                    outflow = outflow + sum(Plink{rcls,scls}(inode,:));
                end
            end
            if outflow < GlobalConstants.FineTol
                model.nodes{inode}.setRouting(model.classes{rcls}, RoutingStrategy.DISABLED);
            end
        end
    end
end

model.link(P);
% link() installs the probabilistic split; the declared strategy replaces it on
% the dispatch (node, class), whose only arcs are the group's targets
for r = 1:size(routedGroupSites,1)
    model.nodes{routedGroupSites(r,1)}.setRouting(model.classes{routedGroupSites(r,2)}, cgroupStrategy(routedGroupSites(r,3)));
end

% Admission constraint on the server station -- see _kb/06-solver-catalog.md (LN section)
if isfield(lqn,'lincon') && idx <= size(lqn.lincon,1) && ~isempty(lqn.lincon{idx,1})
    Aelem = lqn.lincon{idx,1};
    Alayer = zeros(size(Aelem,1), length(model.classes));
    if ishostlayer
        constrainedIdx = lqn.tasksof{idx}; % column j of Aelem is task constrainedIdx(j)
    else
        constrainedIdx = lqn.entriesof{idx}; % column j of Aelem is entry constrainedIdx(j)
    end
    for j = 1:length(constrainedIdx)
        if ishostlayer
            % a task occupies the host through the classes of its activities
            layerClasses = aidxClass(lqn.actsof{constrainedIdx(j)});
        else
            % an entry is occupied by the classes of the calls that target it
            layerClasses = cidxClass(intersect(find(lqn.callpair(:,2) == constrainedIdx(j))', 1:length(cidxClass)));
        end
        for k = 1:length(layerClasses)
            if ~isempty(layerClasses{k})
                Alayer(:, layerClasses{k}.index) = Alayer(:, layerClasses{k}.index) + Aelem(:,j);
            end
        end
    end
    if any(Alayer(:))
        % One region spanning every replica: the constraint models a passive
        % resource of the server as a whole (a semaphore, a connection pool),
        % so replicas share the tokens rather than each holding a private copy
        fcr = model.addRegion(serverStation(1:nreplicas));
        fcr.setConstraint(Alayer, lqn.lincon{idx,2});
    end
end

% Service-rate dependence on a server station -- see _kb/06-solver-catalog.md (LN section)
if isfield(lqn,'lldscaling') || isfield(lqn,'cdscaling') || isfield(lqn,'jdscaling') || isfield(lqn,'pools')
    R = length(model.classes);
    for sidx = idxSet(:)'
        hasld = isfield(lqn,'lldscaling') && sidx <= numel(lqn.lldscaling) && ~isempty(lqn.lldscaling{sidx});
        hascd = isfield(lqn,'cdscaling') && sidx <= numel(lqn.cdscaling) && ~isempty(lqn.cdscaling{sidx});
        hasjd = isfield(lqn,'jdscaling') && sidx <= numel(lqn.jdscaling) && ~isempty(lqn.jdscaling{sidx});
        haspools = isfield(lqn,'pools') && sidx <= numel(lqn.pools) && ~isempty(lqn.pools{sidx});
        if ~(hasld || hascd || hasjd || haspools)
            continue
        end
        sishost = sidx <= lqn.nhosts;
        if sishost
            operandIdx = lqn.tasksof{sidx};
        else
            operandIdx = lqn.entriesof{sidx};
        end
        cols = cell(1,length(operandIdx));
        for j = 1:length(operandIdx)
            if sishost
                % a task occupies the host through the classes of its activities
                layerClasses = aidxClass(lqn.actsof{operandIdx(j)});
            else
                % an entry is occupied by the classes of the calls that target it
                layerClasses = cidxClass(intersect(find(lqn.callpair(:,2) == operandIdx(j))', 1:length(cidxClass)));
            end
            for k = 1:length(layerClasses)
                if ~isempty(layerClasses{k})
                    cols{j}(end+1) = layerClasses{k}.index;
                end
            end
        end
        for m = 1:nreplicas
            if hasld
                srvStation{sidx}{m}.setLimitedLoadDependence(lqn.lldscaling{sidx});
            end
            if hascd
                % beta_{i,r} is product-form only while an operand maps to a single
                % class; where it aggregates several, the same scaling is emitted as
                % a joint dependence, which is numerically identical but not exact
                cdh = lqn_dep_layer_handle(lqn.cdscaling{sidx}, cols, R, model);
                cdpeak = layerPeak(lqn.cdscalingpeak{sidx}, cols, R);
                if all(cellfun(@(c) numel(c) <= 1, cols))
                    srvStation{sidx}{m}.setLimitedClassDependence(cdh, cdpeak);
                else
                    srvStation{sidx}{m}.setLimitedJointDependence(cdh, cdpeak);
                end
            end
            if hasjd
                srvStation{sidx}{m}.setLimitedJointDependence(lqn_dep_layer_handle(lqn.jdscaling{sidx}, cols, R, model), layerPeak(lqn.jdscalingpeak{sidx}, cols, R));
            end
            if haspools
                % A compatibility declaration IS a rate law: the pools clear
                % mu(n) of SN_COMPAT_RATE, order independent at every integer
                % state, which is the REDUNDANCY reading -- every server
                % compatible with a present job works on it and the first to
                % finish cancels the rest. SN_COMPAT_SCALING divides mu(n) by
                % peak*min(1,N/S), which is exactly what the solver's own
                % multiserver term contributes, so the station clears mu(n) and
                % eta is above one at low occupancy. The lowering is to a JOINT
                % dependence, hence an approximation in the layer: see
                % _kb/06-solver-catalog.md (LN section) for why the exact OI
                % analyzer cannot serve a class-switching layer.
                %
                % The PEAK is the rate the pools clear with every server active,
                % sum_t counts*rates: a rate-scaled station reports U = T*S/peak,
                % so leaving it at 1 drops the division by the multiplicity and
                % reports a fully-compatible pool at S times the utilization of
                % the plain multiserver it should reproduce.
                pl = lqn.pools{sidx};
                etaPool = @(nop) sn_compat_scaling(pl.compat, pl.counts, pl.rates, nop);
                poolPeak = sn_compat_peak(pl.counts, pl.rates);
                srvStation{sidx}{m}.setLimitedJointDependence(lqn_dep_layer_handle(etaPool, cols, R, model), layerPeak(poolPeak*ones(1,length(cols)), cols, R));
            end
        end
    end
end

self.ensemble{idx} = model;

    function [P, curClass, jobPos, curStations] = recurActGraph(P, tidx_caller, aidx, curClass, jobPos, curStations)
        jobPosKey(aidx) = jobPos;
        curClassKey{aidx} = curClass;
        curStationKey{aidx} = curStations;
        nextaidxs = find(lqn.graph(aidx,:)); % these include the called entries
        if ~isempty(nextaidxs)
            isNextPrecFork(aidx) = any(isPostAndAct(nextaidxs)); % indexed on aidx to avoid losing it during the recursion
        end
        % pre-fork state, captured at the first branch so that calls this
        % activity issues before the fork stay sequential
        forkSaved = false;
        forkSaveCurClass = curClass;
        forkSaveJobPos = jobPos;
        forkSaveStations = curStations;

        for nextaidx = nextaidxs % for all successor activities
            if ~isempty(nextaidx)
                % Restore pre-fork state for each branch iteration
                if isNextPrecFork(aidx)
                    if ~forkSaved
                        if isPostAndAct(nextaidx)
                            forkSaved = true;
                            forkSaveCurClass = curClass;
                            forkSaveJobPos = jobPos;
                            forkSaveStations = curStations;
                        end
                    else
                        curClass = forkSaveCurClass;
                        jobPos = forkSaveJobPos;
                        curStations = forkSaveStations;
                    end
                end
                isLoop = false;
                % in the activity graph, the following if is entered only
                % by an edge that is the return from a LOOP activity
                if (lqn.graph(aidx,nextaidx) ~= lqn.dag(aidx,nextaidx))
                    isLoop = true;
                end
                if ~(lqn.parent(aidx) == lqn.parent(nextaidx)) % if different parent task
                    % if the successor activity is an entry of another task, this is a call
                    cidx = matchrow(lqn.callpair,[aidx,nextaidx]); % find the call index
                    switch lqn.calltype(cidx)
                        case CallType.ASYNC
                            % Async calls don't modify caller routing - caller continues immediately without blocking.
                            % Arrival rate at destination is handled via arvproc_classes_updmap (lines 170-195, 230-248).
                        case CallType.SYNC
                            [P, jobPos, curClass, curStations] = routeSynchCall(P, jobPos, curClass, curStations);
                        case CallType.FWD
                            % see _kb/06-solver-catalog.md (LN section) for rationale
                    end
                else
                    % at this point, we have processed all calls, let us do the
                    % activities local to the task next
                    if isempty(intersect(lqn.eshift+(1:lqn.nentries), nextaidxs))
                        % if next activity is not an entry
                        jobPos = jobPosKey(aidx);
                        curClass = curClassKey{aidx};
                        curStations = curStationKey{aidx};
                    else
                        if ismember(nextaidxs(find(nextaidxs==nextaidx)-1), lqn.eshift+(1:lqn.nentries))
                            curClassC = curClass;
                        end
                        jobPos = atClient;
                        curClass = curClassC;
                        curStations = {};
                    end
                    % station of the processor the next activity runs on, empty if that processor is not a server of this layer
                    hstn = serversFor(lqn.parent(lqn.parent(nextaidx)));
                    if jobPos == atClient % at client node
                        if ~isempty(hstn)
                            if ~iscachelayer
                                for m=1:length(hstn)
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(clientDelay, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, hstn{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];
                                            P{curClass, forkClass}(clientDelay,joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,hstn{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(clientDelay,hstn{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    hstn{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                    % A SetupTask's cold start is NOT wired into the
                                    % layer station any more. Doing so routed the layer
                                    % through the open M/G/1-with-setup QBD, which reads
                                    % the idle period from the Poisson rate 1/X and so
                                    % powers the thread down far more often than a closed
                                    % layer does -- 12.67% below LDES on lqn_setup -- and
                                    % it charged the setup to the ACTIVITY, where it is
                                    % not host demand. It is charged to the entry instead,
                                    % with probability p: see lqn_setup_charge.
                                end
                                jobPos = atServer;
                                curStations = hstn;
                                curClass = aidxClass{nextaidx};
                                self.servt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, model.getNodeIndex(hstn{1}), aidxClass{nextaidx}.index];
                                % see _kb/06-solver-catalog.md (LN section) for rationale
                                if ~isempty(aidxThinkClass{nextaidx})
                                    for m=1:length(hstn)
                                        P{curClass, aidxThinkClass{nextaidx}}(hstn{m}, clientDelay) = 1.0;
                                    end
                                    curClass = aidxThinkClass{nextaidx};
                                    jobPos = atClient;
                                    curStations = {};
                                end
                            else
                                P{curClass, aidxClass{nextaidx}}(clientDelay,cacheNode) = full(lqn.graph(aidx,nextaidx));

                                cacheNode.setReadItemEntry(aidxClass{nextaidx},lqn.itemproc{aidx},lqn.nitems(aidx));
                                lqn.hitmissaidx = find(lqn.graph(nextaidx,:));
                                lqn.hitaidx = lqn.hitmissaidx(1);
                                lqn.missaidx = lqn.hitmissaidx(2);

                                cacheNode.setHitClass(aidxClass{nextaidx},aidxClass{lqn.hitaidx});
                                cacheNode.setMissClass(aidxClass{nextaidx},aidxClass{lqn.missaidx});

                                if hasRetrievalCache
                                    % see _kb/06-solver-catalog.md (LN section) for rationale
                                    retrievalWiring = struct('readClass', aidxClass{nextaidx}, ...
                                        'missClass', aidxClass{lqn.missaidx}, 'missaidx', lqn.missaidx);
                                end

                                jobPos = atCache; % cache
                                curStations = {};
                                curClass = aidxClass{nextaidx};
                                %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, lqn.hitaidx, 3, 3, aidxClass{nextaidx}.index, aidxClass{lqn.hitaidx}.index];
                                %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, lqn.missaidx, 3, 3, aidxClass{nextaidx}.index, aidxClass{lqn.missaidx}.index];
                            end
                        else % the processor of the next activity is not a server of this layer
                            if isNextPrecFork(aidx)
                                % if next activity is a post-and
                                P{curClass, curClass}(clientDelay, forkNode) = 1.0;
                                f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                forkClassStack(end+1) = curClass.index;
                                P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, clientDelay) = 1.0;
                            else
                                if isPreAndAct(aidx)
                                    % before entering the job we go back to the entry class at the last fork
                                    forkClass = model.classes{forkClassStack(end)};
                                    forkClassStack(end) = [];
                                    P{curClass, forkClass}(clientDelay,joinNode) = 1.0;
                                    P{forkClass, aidxClass{nextaidx}}(joinNode,clientDelay) = 1.0;
                                else
                                    P{curClass, aidxClass{nextaidx}}(clientDelay,clientDelay) = full(lqn.graph(aidx,nextaidx));
                                end
                            end
                            jobPos = atClient;
                            curStations = {};
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.servtproc{nextaidx});
                            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    elseif jobPos == atServer || jobPos == atCache % at server station
                        if ~isempty(hstn)
                            if jobPos == atCache
                                curClass = aidxClass{nextaidx};
                                for m=1:length(hstn)
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(cacheNode, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, hstn{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];

                                            P{curClass, forkClass}(cacheNode,joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,hstn{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(cacheNode,hstn{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    hstn{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                    %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, nextaidx, 3, 2, aidxClass{nextaidx}.index, aidxClass{nextaidx}.index];
                                end
                            else
                                for m=1:length(hstn)
                                    fromStn = curStations{min(m,numel(curStations))};
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(fromStn, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, hstn{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];
                                            P{curClass, forkClass}(fromStn,joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,hstn{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(fromStn,hstn{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    hstn{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                end
                            end
                            jobPos = atServer;
                            curStations = hstn;
                            curClass = aidxClass{nextaidx};
                            self.servt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, model.getNodeIndex(hstn{1}), aidxClass{nextaidx}.index];
                        else
                            for m=1:numel(curStations)
                                fromStn = curStations{m};
                                if isNextPrecFork(aidx)
                                    % if next activity is a post-and
                                    P{curClass, curClass}(fromStn, forkNode) = 1.0;
                                    f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                    forkClassStack(end+1) = curClass.index;
                                    P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                    P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, clientDelay) = 1.0;
                                else
                                    if isPreAndAct(aidx)
                                        % before entering the job we go back to the entry class at the last fork
                                        forkClass = model.classes{forkClassStack(end)};
                                        forkClassStack(end) = [];
                                        P{curClass, forkClass}(fromStn,joinNode) = 1.0;
                                        P{forkClass, aidxClass{nextaidx}}(joinNode,clientDelay) = 1.0;
                                    else
                                        P{curClass, aidxClass{nextaidx}}(fromStn,clientDelay) = full(lqn.graph(aidx,nextaidx));
                                    end
                                end
                            end
                            jobPos = atClient;
                            curStations = {};
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.servtproc{nextaidx});
                            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    end
                    if aidx ~= nextaidx && ~isLoop
                        %% now recursively build the rest of the routing matrix graph
                        [P, curClass, jobPos, curStations] = recurActGraph(P, tidx_caller, nextaidx, curClass, jobPos, curStations);

                        % At this point curClassRec is the last class in the
                        % recursive branch, which we now close with a reply
                        if jobPos == atClient
                            P{curClass, aidxClass{tidx_caller}}(clientDelay,clientDelay) = 1;
                            if ~strcmp(curClass.name(end-3:end),'.Aux')
                                curClass.completes = true;
                            end
                        else
                            for m=1:numel(curStations)
                                P{curClass, aidxClass{tidx_caller}}(curStations{m},clientDelay) = 1;
                            end
                            if ~strcmp(curClass.name(end-3:end),'.Aux')
                                curClass.completes = true;
                            end
                        end
                    end
                end
            end
        end % nextaidx
    end

    function [P, jobPos, curClass, curStations] = routeSynchCall(P, jobPos, curClass, curStations)
        gid = cgroupOfCidx(cidx);
        if gid > 0
            if groupRouted(gid)
                return % the group is one dispatch, already wired at its first member
            end
            [P, jobPos, curClass, curStations, done] = routeGroupCall(P, gid, jobPos, curClass, curStations);
            if done
                groupRouted(gid) = true;
                return
            end
        end
        tstn = serversFor(lqn.parent(lqn.callpair(cidx,2))); % stations of the called task, empty if it is not a server of this layer
        ntgt = length(tstn);
        switch jobPos
            case atClient
                if ~isempty(tstn)
                    % if a call to an entry of a server in this layer
                    % the branch is chosen by callmean alone: a Bernoulli pass can
                    % carry at most one call, more than one needs the geometric
                    % loop through the Aux class. The NTGT replicas only split
                    % each probability, they never relax that bound
                    if callmean(cidx) < 1
                        P{curClass, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 - callmean(cidx);
                        for m=1:ntgt
                            P{curClass, cidxClass{cidx}}(clientDelay,tstn{m}) = callmean(cidx) / ntgt;
                            P{cidxClass{cidx}, cidxClass{cidx}}(tstn{m},clientDelay) = 1.0; % not needed, just to avoid leaving the Aux class disconnected
                        end
                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1.0; % not needed, just to avoid leaving the Aux class disconnected
                    elseif callmean(cidx) == 1
                        for m=1:ntgt
                            P{curClass, cidxClass{cidx}}(clientDelay,tstn{m}) = 1 / ntgt;
                            P{cidxClass{cidx}, cidxClass{cidx}}(tstn{m},clientDelay) = 1.0;
                        end
                    else % callmean(cidx) > 1
                        for m=1:ntgt
                            P{curClass, cidxClass{cidx}}(clientDelay,tstn{m}) = 1 / ntgt;
                            P{cidxClass{cidx}, cidxAuxClass{cidx}}(tstn{m},clientDelay) = 1.0  ;
                            P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,tstn{m}) = (1.0 - 1.0 / callmean(cidx)) / ntgt;
                        end
                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1.0 / (callmean(cidx));
                    end
                    jobPos = atClient;
                    curStations = {};
                    clientDelay.setService(cidxClass{cidx}, Immediate.getInstance());
                    for m=1:ntgt
                        tstn{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(tstn{m}), cidxClass{cidx}.index];
                    end
                    curClass = cidxClass{cidx};
                else
                    % if it is not a call to an entry of a server in this layer
                    if callmean(cidx) < 1
                        % see _kb/06-solver-catalog.md (LN section) for rationale
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1;
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxAuxClass{cidx};
                    elseif callmean(cidx) == 1
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > 1
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1; % the mean number of calls is now embedded in the demand
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;% / (callmean(cidx)/nreplicas); % the mean number of calls is now embedded in the demand
                        curClass = cidxAuxClass{cidx};
                    end
                    jobPos = atClient;
                    curStations = {};
                    clientDelay.setService(cidxClass{cidx}, callservtproc{cidx});
                    self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                end
            case atServer % job at server
                if ~isempty(tstn)
                    % if it is a call to an entry of a server in this layer
                    if callmean(cidx) < 1
                        % The call is skipped with probability 1-callmean; the
                        % two flows merge back in the call class at the client,
                        % as in the atClient branch. Routing the skip into the
                        % call class instead would leave the Aux class with no
                        % inbound arc and its chain without a reference class.
                        for m=1:ntgt
                            fromStn = curStations{min(m,numel(curStations))};
                            P{curClass, cidxAuxClass{cidx}}(fromStn,clientDelay) = 1 - callmean(cidx);
                            P{curClass, cidxClass{cidx}}(fromStn,tstn{m}) = callmean(cidx) / ntgt;
                            P{cidxClass{cidx}, cidxClass{cidx}}(tstn{m},clientDelay) = 1.0;
                            tstn{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        end
                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1.0;
                        % The reply transits the client in the call class, and
                        % sn_refresh_visits drops any (station, class) state whose
                        % rate is NaN, so the class must be declared there
                        clientDelay.setService(cidxClass{cidx}, Immediate.getInstance());
                        jobPos = atClient;
                        curStations = {};
                        curClass = cidxClass{cidx};
                    elseif callmean(cidx) == 1
                        for m=1:ntgt
                            fromStn = curStations{min(m,numel(curStations))};
                            P{curClass, cidxClass{cidx}}(fromStn,tstn{m}) = 1;
                        end
                        if flat
                            % the reply returns the job to the client, which is
                            % where the successor restoration expects it
                            for m=1:ntgt
                                P{cidxClass{cidx}, cidxClass{cidx}}(tstn{m},clientDelay) = 1;
                            end
                            clientDelay.setService(cidxClass{cidx}, Immediate.getInstance());
                            jobPos = atClient;
                            curStations = {};
                        else
                            jobPos = atServer;
                            curStations = tstn;
                        end
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > 1
                        for m=1:ntgt
                            fromStn = curStations{min(m,numel(curStations))};
                            P{curClass, cidxClass{cidx}}(fromStn,tstn{m}) = 1;
                        end
                        if flat
                            % the geometric repeat transits the client between
                            % visits; a self-loop would merge them into one
                            for m=1:ntgt
                                P{cidxClass{cidx}, cidxAuxClass{cidx}}(tstn{m},clientDelay) = 1;
                                P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,tstn{m}) = (1 - 1 / (callmean(cidx))) / ntgt;
                            end
                            P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1 / (callmean(cidx));
                            clientDelay.setService(cidxClass{cidx}, Immediate.getInstance());
                            curClass = cidxClass{cidx};
                        else
                            for m=1:ntgt
                                P{cidxClass{cidx}, cidxClass{cidx}}(tstn{m},tstn{m}) = 1 - 1 / (callmean(cidx));
                                P{cidxClass{cidx}, cidxAuxClass{cidx}}(tstn{m},clientDelay) = 1 / (callmean(cidx));
                            end
                            curClass = cidxAuxClass{cidx};
                        end
                        jobPos = atClient;
                        curStations = {};
                    end
                    for m=1:ntgt
                        tstn{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(tstn{m}), cidxClass{cidx}.index];
                    end
                else
                    % if it is not a call to an entry of a server in this layer
                    % callmean not needed since we switched
                    % to ResidT to model service time at client
                    if callmean(cidx) < 1
                        for m=1:numel(curStations)
                            P{curClass, cidxClass{cidx}}(curStations{m},clientDelay) = 1;
                        end
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxAuxClass{cidx};
                    elseif callmean(cidx) == 1
                        for m=1:numel(curStations)
                            P{curClass, cidxClass{cidx}}(curStations{m},clientDelay) = 1;
                        end
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > 1
                        for m=1:numel(curStations)
                            P{curClass, cidxClass{cidx}}(curStations{m},clientDelay) = 1;
                        end
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxAuxClass{cidx};
                    end
                    jobPos = atClient;
                    curStations = {};
                    clientDelay.setService(cidxClass{cidx}, callservtproc{cidx});
                    self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                end
        end

    end

    function e = entriesOfSet()
        % Entries of every server element of this layer
        e = [];
        for s_ = idxSet(:)'
            e = [e, lqn.entriesof{s_}]; %#ok<AGROW>
        end
    end

    function tf = hostIsServer(tidx_)
        % True when the processor of task TIDX_ is a server of this layer
        tf = ~isempty(serversFor(lqn.parent(tidx_)));
    end

    function st = serversFor(elemIdx)
        % Stations of ELEMIDX when it is a server of this layer, empty otherwise
        if elemIdx >= 1 && elemIdx <= lqn.nidx && ~isempty(srvStation{elemIdx})
            st = srvStation{elemIdx};
        else
            st = {};
        end
    end

    function setAllServers(jobclass, dist)
        % Declare a class at every server station of this layer
        for s_ = idxSet(:)'
            for m_ = 1:length(srvStation{s_})
                srvStation{s_}{m_}.setService(jobclass, dist);
            end
        end
    end

    function [P, jobPos, curClass, curStations, done] = routeGroupCall(P, gid, jobPos, curClass, curStations)
        % A routed group is ONE hop with n destinations. The job is switched into
        % the dispatch class while still at the client, so (client, dispatch)
        % carries exactly the n arcs the strategy chooses among; the 1/n split
        % laid down here is the probabilistic reading a solver without
        % state-dependent routing would see, and is replaced after link().
        done = false;
        tgtStn = {};
        tgtCidx = [];
        for mcidx_ = cgroupMembers{gid}
            mstn_ = serversFor(lqn.parent(lqn.callpair(mcidx_,2)));
            if ~isempty(mstn_)
                tgtStn{end+1} = mstn_{1}; %#ok<AGROW>
                tgtCidx(end+1) = mcidx_; %#ok<AGROW>
            end
        end
        if numel(tgtStn) < 2
            return % not enough of the group is co-resident here to route it
        end
        if jobPos == atClient
            fromStn = clientDelay;
        else
            fromStn = curStations{1};
        end
        dispCls = grpDispatchClass{gid};
        shrCls = grpCallClass{gid};
        P{curClass, dispCls}(fromStn, grpRouter{gid}) = 1.0;
        share_ = 1.0 / numel(tgtStn);
        for m_ = 1:numel(tgtStn)
            P{dispCls, dispCls}(grpRouter{gid}, tgtStn{m_}) = share_;
            P{dispCls, shrCls}(tgtStn{m_}, clientDelay) = 1.0;
            tgtStn{m_}.setService(dispCls, callservtproc{tgtCidx(m_)});
            self.call_classes_updmap{idx}(end+1,:) = [idx, tgtCidx(m_), model.getNodeIndex(tgtStn{m_}), dispCls.index];
        end
        curClass = shrCls;
        jobPos = atClient;
        curStations = {};
        done = true;
    end

    function seedCallService(cidx_)
        % Initial service time of a call class, refined at every iteration
        % through call_classes_updmap. NOTE: under non-flat layering this also
        % seeds the layer's own server for calls whose target is elsewhere,
        % leaving a rate on a class that never visits the station (zero
        % visits). Do not read per-class rates at a station as "served here"
        % without visit-weighting; see sn_has_product_form_not_het_fcfs.
        % Removing the phantom seeding shifts LN fixed points by ~0.5%
        % (mechanism unadjudicated), so it is kept for now.
        if flat
            seedIdx = lqn.parent(lqn.callpair(cidx_,2)); % the called task
        else
            seedIdx = idx;
        end
        minRespT = 0;
        for tidx_act = lqn.actsof{seedIdx}
            minRespT = minRespT + lqn.hostdem{tidx_act}.getMean; % upper bound, uses all activities not just the ones reachable by this entry
        end
        stn_ = serversFor(seedIdx);
        for m_ = 1:length(stn_)
            stn_{m_}.setService(cidxClass{cidx_}, Exp.fitMean(minRespT));
        end
    end

    function st = entryOpenStation(eidx_)
        % Station an open arrival at entry EIDX_ enters, the processor of its
        % task under host layering and the task itself under flat layering
        if flat
            st = serversFor(lqn.parent(eidx_));
        else
            st = serverStation;
        end
    end
end

function peak = layerPeak(peakPerOperand, cols, R)
% PEAK = LAYERPEAK(PEAKPEROPERAND, COLS, R) spread a per-operand peak rate onto layer classes
peak = ones(1,R);
for j = 1:length(cols)
    if ~isempty(cols{j})
        peak(cols{j}) = peakPerOperand(min(j,numel(peakPerOperand)));
    end
end
end

function [ofCidx, members, strategy] = callGroupsByCidx(lqn)
% [OFCIDX, MEMBERS, STRATEGY] = CALLGROUPSBYCIDX(LQN) resolves lqn.callgroups
% from target entries to call indices. OFCIDX(cidx) is the group index of that
% call, 0 when it is dispatched on its own; MEMBERS{g} lists the call indices of
% group g in declaration order and STRATEGY(g) is its RoutingStrategy. A group
% that does not resolve to at least two calls is dropped, so a stale group
% cannot silently rewrite a single call's routing.
ofCidx = zeros(1, lqn.ncalls);
members = {};
strategy = [];
if ~isfield(lqn, 'callgroups') || isempty(lqn.callgroups)
    return
end
for g = 1:numel(lqn.callgroups)
    grp = lqn.callgroups{g};
    found = [];
    for cidx = lqn.callsof{grp.caller}
        if any(lqn.callpair(cidx,2) == grp.targets)
            found(end+1) = cidx; %#ok<AGROW>
        end
    end
    if numel(found) >= 2
        members{end+1} = found; %#ok<AGROW>
        strategy(end+1) = grp.strategy; %#ok<AGROW>
        ofCidx(found) = numel(members);
    end
end
end
