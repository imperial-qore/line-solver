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

% Reference-path interlocking. Under 'refpath' this layer's client is not one
% population per caller but ONE chain running from the reference task down
% through every caller it reaches, so two callers that are the same REF
% customers arriving by two routes stop being counted as two. The tables below
% say which callers are merged into the chain (RPISSPLICED), which synchronous
% calls stop being a surrogate delay at the client and become an explicit
% descent (RPEXPAND), and which entries above or between the callers become hop
% stages on the client Delay (RPISHOP).
% See lqn_ref_routes and _kb/06-solver-catalog.md (LN section).
ncallsz = max(1, lqn.ncalls);
rpOn = false;
rpIsMember = false(lqn.nidx,1);   % task: a caller merged into a chain
rpIsSpliced = false(lqn.nidx,1);  % task: merged AND not the chain head
rpIsHop = false(lqn.nidx,1);      % entry: a stage of the path, no subgraph here
rpExpand = false(ncallsz,1);      % call: re-expanded instead of charged as a delay
rpChildTask = zeros(ncallsz,1);   % call -> callee task when the callee is a member
rpChildEntry = zeros(ncallsz,1);  % call -> callee entry when the callee is a hop
rpRetShare = zeros(ncallsz,1);    % call -> share of the callee's return it takes
rpCallMean = zeros(ncallsz,1);    % call -> mean invocations per execution of its activity
rpHopChildren = cell(lqn.nidx,1); % entry -> its re-expanded calls, descent order
rpHopSeed = zeros(lqn.nidx,1);    % entry -> build-time stage mean
rpMemberInbound = cell(lqn.nidx,1); % task -> the calls entering it
rpHopInbound = cell(lqn.nidx,1);  % entry -> the calls entering it
rpGroupRoots = {};                % per group: the reference task's entries
rpGroupRefTask = [];              % per group: the reference task
rpGroupHeadIsCaller = [];         % per group: the REF task is itself a caller
% Classes minted for the path; the gate and aux of a caller-issued descent are
% the call classes that already exist, so no class is left without a route.
rpRefStage = {};
rpHopCls = cell(lqn.nidx,1);
rpHopRet = cell(lqn.nidx,1);
rpRetCls = cell(lqn.nidx,1);
rpGate = cell(ncallsz,1);
rpAux = cell(ncallsz,1);
rpResume = cell(ncallsz,1);
resolveRefPath();

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
self.refpath_classes_updmap{idx} = zeros(0,7); % [modelidx, nodeidx, classidx, kind, elemidx, sign, fromentry] % reference-path stage means

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
        % A caller merged into a reference chain brings NO population of its
        % own: the customers are the reference task's, and counting this pool as
        % well is the double count refpath exists to remove. Only the CLASS is
        % zeroed -- self.njobs still records the pool, because updateThinkTimes
        % reads max(self.njobs(tidx,:)) across layers and this caller still
        % needs its surrogate delay in every layer that did not merge it.
        rpSplicedHere = rpOn && rpIsSpliced(tidx_caller);
        if rpSplicedHere
            njobsCls = 0;
        else
            njobsCls = njobs;
        end
        aidxClass{tidx_caller} = ClosedClass(model, caller_name, njobsCls, clientDelay);
        clientDelay.setService(aidxClass{tidx_caller}, Disabled.getInstance());
        setAllServers(aidxClass{tidx_caller}, Disabled.getInstance());
        aidxClass{tidx_caller}.completes = false;
        % renormalize residence times using the visits to the task. refreshStruct
        % allows exactly ONE reference class per chain, so a merged caller must
        % yield it to the chain head; ln_layer_refcell then reads the per-caller
        % rate through model.attribute.normclass instead.
        aidxClass{tidx_caller}.setReferenceClass(~rpSplicedHere);
        aidxClass{tidx_caller}.attribute = [LayeredNetworkElement.TASK, tidx_caller];
        model.attribute.tasks(end+1,:) = [aidxClass{tidx_caller}.index, tidx_caller];
        if lqn.isref(tidx_caller)
            clientDelay.setService(aidxClass{tidx_caller}, self.thinkproc{tidx_caller});
        elseif rpSplicedHere
            % the path above this caller is explicit now, so there is no idle
            % time left to fit and no update-map row to refresh it with
            clientDelay.setService(aidxClass{tidx_caller}, Immediate.getInstance());
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

% Reference-path classes. They are minted here, with the caller classes and
% before initRoutingMatrix sizes the routing matrix to the class count: a class
% added afterwards would need P rebuilt around it, as the cache sublayer does.
if rpOn
    buildRefPathClasses();
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

% a priority-scheduling processor orders its callers by task priority
ln_host_task_priorities(lqn, model, idxSet);

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

% Reference-path arcs that no activity graph of this layer could have written:
% the reference stage, the hop stages, and every callee's return split.
if rpOn
    P = routeRefPath(P);
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

% Per-class normalisation task, read by ln_layer_refcell. A merged chain holds
% more than one caller and refreshStruct allows exactly ONE reference class per
% chain, so the invocation rate a residence has to be divided by cannot be found
% through sn.refclass any more; without this every merged caller's residence
% comes out scaled by the product of the call multiplicities along the path.
if rpOn
    model.attribute.normclass = refPathNormClass();
    refPathAudit(P);
end

model.link(P);
% link() installs the probabilistic split; the declared strategy replaces it on
% the dispatch (node, class), whose only arcs are the group's targets
for r = 1:size(routedGroupSites,1)
    model.nodes{routedGroupSites(r,1)}.setRouting(model.classes{routedGroupSites(r,2)}, cgroupStrategy(routedGroupSites(r,3)));
end

% Admission constraint on the server station -- see _kb/06-solver-catalog.md (LN section)
Alayer = zeros(0, length(model.classes));
blayer = [];
if isfield(lqn,'lincon') && idx <= size(lqn.lincon,1) && ~isempty(lqn.lincon{idx,1})
    Aelem = lqn.lincon{idx,1};
    Alayer = zeros(size(Aelem,1), length(model.classes));
    blayer = lqn.lincon{idx,2};
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
end
% A caller merged into a reference chain no longer carries its thread pool as a
% population, so the pool is declared as an admission bound on the layer's
% server instead -- BUT ONLY WHEN THE LAYER SOLVER CAN HONOUR ONE.
% TRAP: sn.regionlincon is read by the CTMC, fluid, SSA and NC loss-network
% analyzers and by NONE of the MVA ones, and LN's default layer solver is
% SolverMVA. That does not make the constraint harmlessly inert there:
% assertLayerSolverSupportsModel REFUSES a layer carrying a region that the
% layer solver's feature set does not cover, so emitting it unconditionally
% made 'refpath' fail to solve at all under the default factory. The pool has a
% second carrier that needs no region -- the thread-acquisition wait
% updateRefPathStages charges at every descent gate and into every hop stage --
% and that is what carries it on the MVA path. There is no double count under a
% region-aware solver either: the customers really are held out there, so the
% wait measured at the task's own layer already reflects the constrained queue.
% See _kb/06-solver-catalog.md (LN section).
if rpOn && self.layerSolverDeclaresRegion()
    [Arp, brp] = refPathPoolRows();
    if ~isempty(Arp)
        Alayer = [Alayer; Arp];
        blayer = [blayer(:); brp(:)];
    end
elseif rpOn
    line_debug(['LN: refpath layer of element %d carries the thread pools of ', ...
        'its merged callers as the measured thread-acquisition wait alone, ', ...
        'because the layer solver declares no Region\n'], idx);
end
if ~isempty(Alayer) && any(Alayer(:))
    % One region spanning every replica: the constraint models a passive
    % resource of the server as a whole (a semaphore, a connection pool),
    % so replicas share the tokens rather than each holding a private copy
    fcr = model.addRegion(serverStation(1:nreplicas));
    fcr.setConstraint(Alayer, blayer);
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
                        % recursive branch, which we now close with a reply.
                        % A caller merged into a reference chain replies to its
                        % OWN return class instead, which hands the customer
                        % back to whichever descent invoked it.
                        replyCls = aidxClass{tidx_caller};
                        if rpOn && rpIsSpliced(tidx_caller)
                            replyCls = rpRetCls{tidx_caller};
                        end
                        if jobPos == atClient
                            P{curClass, replyCls}(clientDelay,clientDelay) = 1;
                            if ~strcmp(curClass.name(end-3:end),'.Aux')
                                curClass.completes = true;
                            end
                        else
                            for m=1:numel(curStations)
                                P{curClass, replyCls}(curStations{m},clientDelay) = 1;
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
        if rpOn && rpExpand(cidx)
            % The callee is on the reference path into this layer: its time is
            % not a surrogate delay on the client any more, the customer
            % descends into it. resolveRefPath has already refused the layer if
            % such a callee were served by one of this layer's own stations.
            [P, jobPos, curClass, curStations] = routeRefPathDescent(P, cidx, jobPos, curClass, curStations);
            return
        end
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

    function resolveRefPath()
        % Decide whether this layer's client chains are merged along the
        % reference path, and fill the tables the three refpath passes read.
        % Everything here is a pure query over LQN plus the layer's own server
        % set: no class exists yet when it runs, because the first pass needs
        % the answer to size the caller populations.
        if ~strcmp(self.interlockMethod, 'refpath') || flat || iscachelayer || hasfork || hasjoin
            return
        end
        % A Fork in the layer flattens every positive alpha in
        % sn_refresh_visits, which would destroy the geometric-loop visit ratios
        % the descent is built from; a replicated element carries a per-replica
        % njobs convention that does not compose with a shared chain population.
        if reduceFanout || nreplicas > 1
            return
        end
        for t_ = callers(:)'
            if lqn.repl(t_) > 1 || any(self.singleReplicaTasks == t_)
                return
            end
        end
        maxpaths = 32;
        if isfield(self.options.config,'interlock_maxpaths') && ~isempty(self.options.config.interlock_maxpaths)
            maxpaths = self.options.config.interlock_maxpaths;
        end
        scope = 'merging';
        if isfield(self.options.config,'interlock_refpath_scope') && ~isempty(self.options.config.interlock_refpath_scope)
            scope = lower(self.options.config.interlock_refpath_scope);
        end
        [R, why] = lqn_ref_routes(lqn, callers, maxpaths, idxSet);
        if ~isempty(why)
            % Named, because the discontinuity is otherwise invisible: a
            % one-line model edit can cross the route bound or add a second
            % reference task and every number in the layer moves.
            line_warning(mfilename, ['layer ''%s'' falls back to the ''ilrate'' interlock ' ...
                'method: %s.\n'], lqn.hashnames{idx}, why);
            return
        end
        keep = [];
        for g_ = 1:numel(R)
            wanted = numel(R(g_).members) >= 2;
            if ~wanted && strcmp(scope,'all') && ~R(g_).headIsCaller && ~isempty(R(g_).members)
                % a lone caller below a reference task: no merge, but the path
                % above it becomes explicit instead of a fitted idle time
                wanted = true;
            end
            if ~wanted
                continue
            end
            % The merge only pays when the chain population it installs is
            % STRICTLY below the caller pools it replaces. mult(REF) is uncapped
            % by design, so a reference task carrying as many threads as its
            % callers hold between them leaves the layer with exactly the client
            % population it had -- while the rate discount of 'ilrate' has been
            % switched off along with the double count it was correcting, which
            % makes refpath strictly worse than doing nothing there. Measured on
            % 17-external-interlock, mult(REF) = 2 against two callers of one
            % thread each: the correction goes to zero and the discount with it.
            % The group is handed back to 'ilrate' rather than merged for free.
            % Dropping ONE group is safe where refusing one caller is not:
            % lqn_ref_routes has already rejected the layer if any caller
            % descends from two reference tasks, so the groups partition the
            % callers and a dropped group keeps its ordinary per-caller chains.
            pool = 0;
            for t_ = R(g_).members(:)'
                pool = pool + refPathCallerPool(t_);
            end
            chainPop = mult(R(g_).reftask) * lqn.repl(R(g_).reftask);
            if ~(pool > chainPop)
                line_debug(['LN refpath: layer %s keeps a caller group on ''ilrate'', ', ...
                    'reference task %s holds %g thread(s) against %g in the %d caller(s) ', ...
                    'it would merge, so the merge would remove no customer\n'], ...
                    lqn.hashnames{idx}, lqn.hashnames{R(g_).reftask}, ...
                    chainPop, pool, numel(R(g_).members));
                continue
            end
            keep(end+1) = g_; %#ok<AGROW>
        end
        if isempty(keep)
            return
        end
        % Stage the whole decision, then commit: a refusal found on the second
        % group must not leave the first one half-applied.
        sIsMember = false(lqn.nidx,1);
        sIsSpliced = false(lqn.nidx,1);
        sIsHop = false(lqn.nidx,1);
        sExpand = false(ncallsz,1);
        sChildTask = zeros(ncallsz,1);
        sChildEntry = zeros(ncallsz,1);
        sRetShare = zeros(ncallsz,1);
        sCallMean = zeros(ncallsz,1);
        sHopChildren = cell(lqn.nidx,1);
        sHopSeed = zeros(lqn.nidx,1);
        sMemberInbound = cell(lqn.nidx,1);
        sHopInbound = cell(lqn.nidx,1);
        sRoots = {};
        sRefTask = [];
        sHeadIsCaller = [];
        for g_ = keep
            G = R(g_);
            nE = numel(G.entries);
            % Does this DAG node lead to a caller of the layer? Reverse
            % topological order, so every child is already decided.
            leads = G.ismember(:)';
            for i_ = nE:-1:1
                for k_ = find(G.calls(:,2)' == i_)
                    if leads(G.calls(k_,3))
                        leads(i_) = true;
                    end
                end
            end
            for t_ = G.members(:)'
                if ~callerHasClass(t_)
                    return % the caller has no client class to enter
                end
            end
            if G.headIsCaller && ~callerHasClass(G.reftask)
                return
            end
            if ~G.headIsCaller && ~isfinite(lqn.maxmult(G.reftask))
                return % an infinite reference pool has no chain population
            end
            for i_ = 1:nE
                if ~leads(i_) || G.ismember(i_)
                    continue
                end
                % A hop whose task is a server of this layer would charge at the
                % client the work the layer's own station exists to serve, and
                % the customer would be in two places at once.
                if any(idxSet == G.etask(i_)) || lqn.repl(G.etask(i_)) > 1
                    return
                end
            end
            for k_ = 1:size(G.calls,1)
                if ~leads(G.calls(k_,3))
                    continue
                end
                if any(idxSet == G.etask(G.calls(k_,3)))
                    return % the callee is served by this layer's own station
                end
            end
            % Commit this group into the staging tables
            head = G.reftask;
            sIsMember(G.members) = true;
            for t_ = G.members(:)'
                sIsSpliced(t_) = ~(G.headIsCaller && t_ == head);
            end
            inTotMember = zeros(lqn.nidx,1);
            inTotHop = zeros(lqn.nidx,1);
            for k_ = 1:size(G.calls,1)
                toPos = G.calls(k_,3);
                if ~leads(toPos)
                    continue
                end
                cidx_ = G.calls(k_,1);
                fromEidx = G.entries(G.calls(k_,2));
                sExpand(cidx_) = true;
                sCallMean(cidx_) = lqn.callproc{cidx_}.getMean;
                sHopChildren{fromEidx}(end+1) = cidx_;
                if G.ismember(toPos)
                    tv = G.etask(toPos);
                    sChildTask(cidx_) = tv;
                    sMemberInbound{tv}(end+1) = cidx_;
                    inTotMember(tv) = inTotMember(tv) + G.calls(k_,5);
                else
                    ev = G.entries(toPos);
                    sChildEntry(cidx_) = ev;
                    sHopInbound{ev}(end+1) = cidx_;
                    inTotHop(ev) = inTotHop(ev) + G.calls(k_,5);
                end
                sRetShare(cidx_) = G.calls(k_,5); % normalised below
            end
            for k_ = 1:size(G.calls,1)
                cidx_ = G.calls(k_,1);
                if ~sExpand(cidx_)
                    continue
                end
                if sChildTask(cidx_) > 0
                    tot = inTotMember(sChildTask(cidx_));
                else
                    tot = inTotHop(sChildEntry(cidx_));
                end
                if tot > GlobalConstants.FineTol
                    sRetShare(cidx_) = sRetShare(cidx_) / tot;
                else
                    sRetShare(cidx_) = 0;
                end
            end
            for i_ = 1:nE
                if leads(i_) && ~G.ismember(i_)
                    ev = G.entries(i_);
                    sIsHop(ev) = true;
                    sHopSeed(ev) = hostDemandOfEntry(ev);
                end
            end
            if ~G.headIsCaller
                % EVERY entry of the reference task becomes a stage, not only
                % those that reach this layer. The reference task splits its
                % cycle across all of them -- the same 1/ncaller_entries the
                % second pass lays down for a caller -- so leaving the others
                % out would shrink the denominator and make every cycle reach
                % this layer, over-driving the chain by exactly the ratio of
                % leading entries to entries. An entry that leads nowhere has no
                % descents to re-expand, so its stage is simply its own
                % residence and it returns straight to the reference stage.
                for ev = lqn.entriesof{G.reftask}(:)'
                    if ~sIsHop(ev)
                        sIsHop(ev) = true;
                        sHopSeed(ev) = hostDemandOfEntry(ev);
                    end
                end
            end
            sRoots{end+1} = lqn.entriesof{head}(:)'; %#ok<AGROW>
            sRefTask(end+1) = head; %#ok<AGROW>
            sHeadIsCaller(end+1) = G.headIsCaller; %#ok<AGROW>
        end
        rpIsMember = sIsMember;
        rpIsSpliced = sIsSpliced;
        rpIsHop = sIsHop;
        rpExpand = sExpand;
        rpChildTask = sChildTask;
        rpChildEntry = sChildEntry;
        rpRetShare = sRetShare;
        rpCallMean = sCallMean;
        rpHopChildren = sHopChildren;
        rpHopSeed = sHopSeed;
        rpMemberInbound = sMemberInbound;
        rpHopInbound = sHopInbound;
        rpGroupRoots = sRoots;
        rpGroupRefTask = sRefTask;
        rpGroupHeadIsCaller = sHeadIsCaller;
        rpOn = true;
        line_debug('LN refpath: layer %s merges %d caller(s) into %d reference chain(s)', ...
            lqn.hashnames{idx}, sum(rpIsMember), numel(rpGroupRefTask));
    end

    function [P, jobPos, curClass, curStations] = routeRefPathDescent(P, cidx_, jobPos, curClass, curStations)
        % One explicit descent into the callee of CIDX_, replacing the surrogate
        % Delay that would otherwise charge the callee's whole residence at the
        % client. The multiplicity is a geometric loop per edge, never a single
        % flat loop of the product: with t0 --x2--> t1 --x3--> c, t1 runs twice
        % and c six times per reference cycle, and a flattened loop visits t1
        % once. The m < 1 branch is not optional either -- 1 - 1/m goes negative
        % below one, and model.setChecks(false) would let that reach dtmc_solve.
        gate = rpGate{cidx_};
        res = rpResume{cidx_};
        aux = rpAux{cidx_};
        m_ = rpCallMean(cidx_);
        if rpChildTask(cidx_) > 0
            enterCls = aidxClass{rpChildTask(cidx_)};
        else
            enterCls = rpHopCls{rpChildEntry(cidx_)};
        end
        if jobPos == atClient
            P{curClass, gate}(clientDelay, clientDelay) = 1;
        else
            for s_ = 1:numel(curStations)
                P{curClass, gate}(curStations{s_}, clientDelay) = 1;
            end
        end
        if m_ == 1
            P{gate, enterCls}(clientDelay, clientDelay) = 1;
            cont = res;
        elseif m_ > 1
            P{gate, enterCls}(clientDelay, clientDelay) = 1;
            P{res, enterCls}(clientDelay, clientDelay) = 1 - 1/m_;
            P{res, aux}(clientDelay, clientDelay) = 1/m_;
            cont = aux;
        else
            P{gate, enterCls}(clientDelay, clientDelay) = m_;
            P{gate, aux}(clientDelay, clientDelay) = 1 - m_;
            P{res, aux}(clientDelay, clientDelay) = 1;
            cont = aux;
        end
        jobPos = atClient;
        curStations = {};
        curClass = cont;
    end

    function P = routeRefPath(P)
        % The part of the reference path that has no activity graph in this
        % layer: the reference stage itself, the hop stages between the
        % reference task and the callers, and the return arcs that hand a
        % callee's reply back to whichever descent invoked it.
        for g_ = 1:numel(rpGroupRefTask)
            if rpGroupHeadIsCaller(g_)
                continue % the head's own activity graph carries the descent
            end
            roots = rpGroupRoots{g_};
            roots = roots(rpIsHop(roots));
            if isempty(roots)
                continue
            end
            % A reference task picks among its own entries, so this is a
            % PROBABILISTIC split, matching the 1/nentries the second pass lays
            % down for a caller class. The descents below a hop are a SERIES:
            % one reference cycle traverses every one of them.
            for r_ = roots
                P{rpRefStage{g_}, rpHopCls{r_}}(clientDelay, clientDelay) = 1/numel(roots);
                P{rpHopRet{r_}, rpRefStage{g_}}(clientDelay, clientDelay) = 1;
            end
        end
        for ev = find(rpIsHop)'
            cur = rpHopCls{ev};
            kids = rpHopChildren{ev};
            for k_ = 1:numel(kids)
                [P, ~, cur, ~] = routeRefPathDescent(P, kids(k_), atClient, cur, {});
            end
            P{cur, rpHopRet{ev}}(clientDelay, clientDelay) = 1;
        end
        % Return arcs, split by the share of the callee's invocations each
        % descent contributes. Flow into the callee from parent i is v(p_i) on
        % the descent arc plus v(u)*w_i*(1-1/m_i) on the resume back-arc, which
        % is v(p_i)*m_i: the shares are exact, not a heuristic.
        for tv = find(rpIsSpliced)'
            for cidx_ = rpMemberInbound{tv}
                if rpRetShare(cidx_) > 0
                    P{rpRetCls{tv}, rpResume{cidx_}}(clientDelay, clientDelay) = rpRetShare(cidx_);
                end
            end
        end
        for ev = find(rpIsHop)'
            for cidx_ = rpHopInbound{ev}
                if rpRetShare(cidx_) > 0
                    P{rpHopRet{ev}, rpResume{cidx_}}(clientDelay, clientDelay) = rpRetShare(cidx_);
                end
            end
        end
    end

    function buildRefPathClasses()
        % Mint the classes the reference path needs. A descent issued from a
        % caller's own activity graph REUSES the call classes that already exist
        % as its gate and its geometric auxiliary: minting fresh ones would
        % leave those two without a single inbound arc, and link() gives a class
        % with no routing row RAND routing at every node.
        for g_ = 1:numel(rpGroupRefTask)
            if rpGroupHeadIsCaller(g_)
                rpRefStage{g_} = []; %#ok<AGROW> % the caller class IS the head
                continue
            end
            rt = rpGroupRefTask(g_);
            njobsRef = mult(rt) * lqn.repl(rt);
            cls = ClosedClass(model, [lqn.hashnames{rt},'.RefPath'], njobsRef, clientDelay);
            clientDelay.setService(cls, self.thinkproc{rt});
            setAllServers(cls, Disabled.getInstance());
            cls.completes = false;
            cls.setReferenceClass(true);
            cls.attribute = [-1, rt];
            rpRefStage{g_} = cls; %#ok<AGROW>
        end
        for ev = find(rpIsHop)'
            cls = ClosedClass(model, [lqn.hashnames{ev},'.RefHop'], 0, clientDelay);
            if rpHopSeed(ev) > GlobalConstants.FineTol
                clientDelay.setService(cls, Exp.fitMean(rpHopSeed(ev)));
            else
                clientDelay.setService(cls, Immediate.getInstance());
            end
            setAllServers(cls, Disabled.getInstance());
            cls.completes = false;
            cls.attribute = [-1, ev];
            rpHopCls{ev} = cls;
            % The stage charges the entry's own residence, MINUS the descents
            % the chain re-expands below it, PLUS the wait its task's threads
            % queue for. updateRefPathStages sums these every iteration.
            self.refpath_classes_updmap{idx}(end+1,:) = [idx, 1, cls.index, 1, ev, 1, 0];
            for c_ = rpHopChildren{ev}
                % Column 7 is the entry the call is issued FROM, recorded here
                % because the build knows it exactly. Inverting actsof instead
                % would be ambiguous: actsof{eidx} is a reachability closure, so
                % two entries of one task that converge on a shared activity
                % both contain it and the inverse has no single answer.
                self.refpath_classes_updmap{idx}(end+1,:) = [idx, 1, cls.index, 2, c_, -1, ev];
            end
            self.refpath_classes_updmap{idx}(end+1,:) = [idx, 1, cls.index, 3, lqn.parent(ev), 1, 0];
            rcls = ClosedClass(model, [lqn.hashnames{ev},'.RefHopRet'], 0, clientDelay);
            clientDelay.setService(rcls, Immediate.getInstance());
            setAllServers(rcls, Disabled.getInstance());
            rcls.completes = false;
            rcls.attribute = [-1, ev];
            rpHopRet{ev} = rcls;
        end
        for tv = find(rpIsSpliced)'
            cls = ClosedClass(model, [lqn.hashnames{tv},'.RefRet'], 0, clientDelay);
            clientDelay.setService(cls, Immediate.getInstance());
            setAllServers(cls, Disabled.getInstance());
            cls.completes = false;
            cls.attribute = [-1, tv];
            rpRetCls{tv} = cls;
        end
        for cidx_ = find(rpExpand)'
            if numel(cidxClass) >= cidx_ && ~isempty(cidxClass{cidx_})
                rpGate{cidx_} = cidxClass{cidx_};
                clientDelay.setService(rpGate{cidx_}, Immediate.getInstance());
                if numel(cidxAuxClass) >= cidx_ && ~isempty(cidxAuxClass{cidx_})
                    rpAux{cidx_} = cidxAuxClass{cidx_};
                end
            else
                cls = ClosedClass(model, [lqn.callhashnames{cidx_},'.RefGate'], 0, clientDelay);
                clientDelay.setService(cls, Immediate.getInstance());
                setAllServers(cls, Disabled.getInstance());
                cls.completes = false;
                cls.attribute = [-1, cidx_];
                rpGate{cidx_} = cls;
                if rpCallMean(cidx_) ~= 1
                    acls = ClosedClass(model, [lqn.callhashnames{cidx_},'.RefGate.Aux'], 0, clientDelay);
                    clientDelay.setService(acls, Immediate.getInstance());
                    setAllServers(acls, Disabled.getInstance());
                    acls.completes = false;
                    acls.attribute = [-1, cidx_];
                    rpAux{cidx_} = acls;
                end
            end
            if rpChildTask(cidx_) > 0
                % A MEMBER is entered through this gate, and a member has no hop
                % stage of its own: its subgraph is explicit in this layer, so
                % every activity residence and every call residence it owns is
                % already charged. The ONE thing that is not is the wait for one
                % of its own threads, which lives in its own task layer and is
                % what its pool means. Without this row a merged caller's pool
                % reaches the layer through the admission region alone, and that
                % region is inert under the default MVA layer solver -- so the
                % pool would silently bound nothing at all. A gate descending
                % into a HOP gets no row: the hop stage charges W for its task.
                self.refpath_classes_updmap{idx}(end+1,:) = ...
                    [idx, 1, rpGate{cidx_}.index, 3, rpChildTask(cidx_), 1, 0];
            end
            scls = ClosedClass(model, [lqn.callhashnames{cidx_},'.RefResume'], 0, clientDelay);
            clientDelay.setService(scls, Immediate.getInstance());
            setAllServers(scls, Disabled.getInstance());
            scls.completes = false;
            scls.attribute = [-1, cidx_];
            rpResume{cidx_} = scls;
        end
    end

    function nc = refPathNormClass()
        % Per class, the TASK class whose invocation rate normalises it. Filled
        % only for the classes of a merged caller: every other class keeps 0 and
        % ln_layer_refcell falls back on the chain's reference class, which is
        % the class it would have picked anyway.
        nc = zeros(1, model.getNumberOfClasses);
        for tv = find(rpIsMember)'
            if isempty(aidxClass{tv})
                continue
            end
            own = aidxClass{tv}.index;
            owned = {aidxClass{tv}};
            for e_ = lqn.entriesof{tv}
                owned{end+1} = aidxClass{e_}; %#ok<AGROW>
            end
            for a_ = lqn.actsof{tv}
                owned{end+1} = aidxClass{a_}; %#ok<AGROW>
                owned{end+1} = aidxThinkClass{a_}; %#ok<AGROW>
                for c_ = lqn.callsof{a_}
                    if numel(cidxClass) >= c_
                        owned{end+1} = cidxClass{c_}; %#ok<AGROW>
                    end
                    if numel(cidxAuxClass) >= c_
                        owned{end+1} = cidxAuxClass{c_}; %#ok<AGROW>
                    end
                end
            end
            for q_ = 1:numel(owned)
                % An open class (an async call, an entry arrival) has no chain
                % reference rate of its own, so it must keep the fallback
                if ~isempty(owned{q_}) && isa(owned{q_}, 'ClosedClass')
                    nc(owned{q_}.index) = own;
                end
            end
        end
    end

    function refPathAudit(P)
        % model.setChecks(false) at the top of this function disables BOTH of
        % link()'s guards, the negative-probability one and the row-sum one, so
        % a mis-wired descent would reach dtmc_solve as garbage with no message
        % at all.
        %
        % ONLY THE CLASSES REFPATH WROTE ARE CHECKED. Sweeping every class would
        % re-audit routing this function has always written, whose partial rows
        % are deliberate -- routeSynchCall leaves several arcs in place purely so
        % an Aux class is not disconnected -- and a flood of warnings about
        % pre-existing structure would bury the one row that is actually wrong.
        if isa(P, 'RoutingMatrix')
            Pc = P.getCell();
        else
            Pc = P;
        end
        K = model.getNumberOfClasses;
        I = length(model.nodes);
        watch = refPathClassIndices();
        for r_ = watch
            if r_ < 1 || r_ > min(K, size(Pc,1))
                continue
            end
            for n_ = 1:I
                tot = 0;
                neg = false;
                for s_ = 1:min(K, size(Pc,2))
                    if isempty(Pc{r_,s_})
                        continue
                    end
                    row = Pc{r_,s_}(n_,:);
                    if any(row < -GlobalConstants.FineTol)
                        neg = true;
                    end
                    tot = tot + sum(row);
                end
                if neg
                    line_warning(mfilename, ['refpath layer ''%s'': class ''%s'' routes ' ...
                        'out of node %d with a negative probability.\n'], ...
                        lqn.hashnames{idx}, model.classes{r_}.name, n_);
                end
                if tot > GlobalConstants.FineTol && abs(tot - 1) > 1e-6
                    line_warning(mfilename, ['refpath layer ''%s'': class ''%s'' leaves ' ...
                        'node %d with total probability %g.\n'], ...
                        lqn.hashnames{idx}, model.classes{r_}.name, n_, tot);
                end
            end
        end
        % Exactly one reference class per chain, which is what refreshStruct
        % requires: two isrefclass members of one chain make refclass(c) a
        % 2-vector and MATLAB throws inside getStruct() with a stack that never
        % mentions LN. Counting them here names the layer instead.
        nref = 0;
        for k_ = 1:numel(model.classes)
            if model.classes{k_}.isrefclass
                nref = nref + 1;
            end
        end
        % One per group whose head is NOT a caller (the minted RefPath class),
        % plus one per caller class that was not spliced into a chain.
        expect = sum(~rpGroupHeadIsCaller);
        for t_ = callers(:)'
            if ~rpIsSpliced(t_) && callerHasClass(t_)
                expect = expect + 1;
            end
        end
        if nref ~= expect
            line_warning(mfilename, ['refpath layer ''%s'' declares %d reference classes ' ...
                'where %d were built; a chain holding two makes refreshStruct throw ' ...
                'inside getStruct().\n'], lqn.hashnames{idx}, nref, expect);
        end
    end

    function ix = refPathClassIndices()
        % Every class the reference path wrote a routing row for, including the
        % call classes a caller-issued descent REUSES as its gate and auxiliary
        % and the caller classes a descent enters.
        ix = zeros(1,0);
        for g_ = 1:numel(rpRefStage)
            if ~isempty(rpRefStage{g_})
                ix(end+1) = rpRefStage{g_}.index; %#ok<AGROW>
            end
        end
        for ev = find(rpIsHop)'
            if ~isempty(rpHopCls{ev}); ix(end+1) = rpHopCls{ev}.index; end %#ok<AGROW>
            if ~isempty(rpHopRet{ev}); ix(end+1) = rpHopRet{ev}.index; end %#ok<AGROW>
        end
        for tv = find(rpIsMember)'
            if ~isempty(rpRetCls{tv}); ix(end+1) = rpRetCls{tv}.index; end %#ok<AGROW>
            if ~isempty(aidxClass{tv}); ix(end+1) = aidxClass{tv}.index; end %#ok<AGROW>
        end
        for c_ = find(rpExpand)'
            if ~isempty(rpGate{c_}); ix(end+1) = rpGate{c_}.index; end %#ok<AGROW>
            if ~isempty(rpAux{c_}); ix(end+1) = rpAux{c_}.index; end %#ok<AGROW>
            if ~isempty(rpResume{c_}); ix(end+1) = rpResume{c_}.index; end %#ok<AGROW>
        end
        ix = unique(ix);
    end

    function [A, b] = refPathPoolRows()
        % One admission row per merged caller: at most maxmult(t) of its threads
        % may hold the layer's server at once. The head needs none, its pool IS
        % the chain population.
        A = zeros(0, length(model.classes));
        b = [];
        for tv = find(rpIsSpliced)'
            cap = full(lqn.maxmult(tv));
            if ~isfinite(cap) || cap <= 0
                continue
            end
            row = zeros(1, length(model.classes));
            if ishostlayer
                % a task occupies the host through the classes of its activities
                for a_ = lqn.actsof{tv}
                    if ~isempty(aidxClass{a_})
                        row(aidxClass{a_}.index) = 1;
                    end
                end
            else
                % a task occupies the server through the calls it makes into it
                srvEntries = entriesOfSet();
                for a_ = lqn.actsof{tv}
                    for c_ = lqn.callsof{a_}
                        if any(lqn.callpair(c_,2) == srvEntries) && numel(cidxClass) >= c_ ...
                                && ~isempty(cidxClass{c_})
                            row(cidxClass{c_}.index) = 1;
                        end
                    end
                end
            end
            if any(row)
                A(end+1,:) = row; %#ok<AGROW>
                b(end+1,1) = cap; %#ok<AGROW>
            end
        end
    end

    function d = hostDemandOfEntry(eidx_)
        % Build-time seed for a hop stage. Iteration 1 runs before any
        % updateLayers, and an Immediate seed leaves the first iterate badly
        % over-driven, so the stage starts at the entry's own host demand.
        d = 0;
        for a_ = lqn.actsof{eidx_}
            if ~isempty(lqn.hostdem{a_})
                m_ = lqn.hostdem{a_}.getMean();
                if isfinite(m_)
                    d = d + m_;
                end
            end
        end
    end

    function tf = callerHasClass(tidx_)
        % The predicate the first pass uses to decide whether caller TIDX_ gets
        % a client class. Mirrored here because the refpath tables are resolved
        % before that pass runs and a member with no class would be unreachable.
        hasDirect = false;
        hasArrival = false;
        isFwdTgt = false;
        if hostIsServer(tidx_)
            if lqn.isref(tidx_)
                hasDirect = true;
            else
                for e_ = lqn.entriesof{tidx_}
                    if any(full(lqn.issynccaller(:, e_))) || any(full(lqn.isasynccaller(:, e_)))
                        hasDirect = true;
                    end
                    if isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && iscell(lqn.arrival) ...
                            && e_ <= length(lqn.arrival) && ~isempty(lqn.arrival{e_})
                        hasArrival = true;
                    end
                    for c_ = 1:lqn.ncalls
                        if full(lqn.calltype(c_)) == CallType.FWD && full(lqn.callpair(c_,2)) == e_
                            isFwdTgt = true;
                            break
                        end
                    end
                end
            end
        end
        tf = (hostIsServer(tidx_) && (hasDirect || hasArrival || isFwdTgt)) ...
            || any(any(lqn.issynccaller(tidx_, entriesOfSet())));
    end

    function n_ = refPathCallerPool(tidx_)
        % The population caller TIDX_ would contribute to this layer WITHOUT the
        % merge. Mirrored from the first pass, including the infinite-server
        % fallbacks, because RESOLVEREFPATH has to compare it against the chain
        % population before any class exists. Kept as a copy rather than
        % factored out of the first pass: that block also WRITES self.njobs and
        % is read by updateThinkTimes, and the two responsibilities should not
        % be entangled for the sake of one comparison.
        if self.njobs(tidx_,idx) ~= 0
            n_ = self.njobs(tidx_,idx);
            return
        end
        if reduceFanout || any(self.singleReplicaTasks == tidx_)
            n_ = mult(tidx_);
        else
            n_ = mult(tidx_) * lqn.repl(tidx_);
        end
        if isinf(n_)
            n_ = sum(mult(find(lqn.taskgraph(:,tidx_)))); %#ok<FNDSB>
            if isinf(n_)
                n_ = min(sum(mult(isfinite(mult)) .* lqn.repl(isfinite(mult))), 1000);
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
