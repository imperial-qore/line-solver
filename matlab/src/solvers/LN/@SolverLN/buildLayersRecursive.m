function  buildLayersRecursive(self, idx, callers, ishostlayer)
lqn = self.lqn;
jobPosKey = zeros(lqn.nidx,1);
curClassKey = cell(lqn.nidx,1);
nreplicas = lqn.repl(idx);
%mult = lqn.mult;
mult = lqn.maxmult; % this removes spare capacity that cannot be used
lqn.mult = mult;
callservtproc = self.callservtproc;
model = Network(lqn.hashnames{idx});
model.setChecks(false); % fast mode
model.attribute = struct('hosts',[],'tasks',[],'entries',[],'activities',[],'calls',[],'serverIdx',0);
if ishostlayer | any(any(lqn.issynccaller(callers, lqn.entriesof{idx}))) | any(any(lqn.isasynccaller(callers, lqn.entriesof{idx}))) %#ok<OR2>
    clientDelay = Delay(model, 'Clients');
    model.attribute.clientIdx = 1;
    model.attribute.serverIdx = 2;
    model.attribute.sourceIdx = NaN;
else
    model.attribute.serverIdx = 1;
    model.attribute.clientIdx = NaN;
    model.attribute.sourceIdx = NaN;
end
serverStation = cell(1,nreplicas);
isfunctionlayer = all(lqn.isfunction(callers)) && ishostlayer;
for m=1:nreplicas
    if m == 1
        serverStation{m} = Queue(model,lqn.hashnames{idx}, lqn.sched(idx));
    else
        serverStation{m} = Queue(model,[lqn.hashnames{idx},'.',num2str(m)], lqn.sched(idx));
    end
    serverStation{m}.setNumberOfServers(mult(idx));
    serverStation{m}.attribute.ishost = ishostlayer;
    serverStation{m}.attribute.idx = idx;
end

iscachelayer = all(lqn.iscache(callers)) && ishostlayer;
if iscachelayer
    cacheNode = Cache(model, lqn.hashnames{callers}, lqn.nitems(callers), lqn.itemcap{callers}, lqn.replacestrat(callers));
end

actsInCaller = lqn.actsof{callers};
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

self.servt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % server classes to update
self.thinkt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % client classes to update
self.actthinkt_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % activity think-time classes to update
self.arvproc_classes_updmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % classes to update in the next iteration for asynch calls
self.call_classes_updmap{idx} = zeros(0,4); % [modelidx, callidx, node, class] % calls classes to update in the next iteration (includes calls in client classes)
self.route_prob_updmap{idx} = zeros(0,7); % [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto] % routing probabilities to update in the next iteration

if ishostlayer
    model.attribute.hosts(end+1,:) = [NaN, model.attribute.serverIdx ];
else
    model.attribute.tasks(end+1,:) = [NaN, model.attribute.serverIdx ];
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
    if ishostlayer
        % Check if the task is a reference task (always create closed class)
        if lqn.isref(tidx_caller)
            hasDirectCallers = true;
        else
            % Check if any entry of this task has sync or async callers
            for eidx = lqn.entriesof{tidx_caller}
                if any(full(lqn.issynccaller(:, eidx))) || any(full(lqn.isasynccaller(:, eidx)))
                    hasDirectCallers = true;
                    break;
                end
                % Also check for open arrivals on this entry
                if isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && ...
                   iscell(lqn.arrival) && eidx <= length(lqn.arrival) && ...
                   ~isempty(lqn.arrival{eidx})
                    hasDirectCallers = true;
                    break;
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
    if (ishostlayer && hasDirectCallers) | any(any(lqn.issynccaller(tidx_caller, lqn.entriesof{idx}))) %#ok<OR2> % if it is only an asynch caller the closed classes are not needed
        if self.njobs(tidx_caller,idx) == 0
            % for each entry of the calling task
            % determine job population
            % this block matches the corresponding calculations in
            % updateThinkTimes
            njobs = mult(tidx_caller)*lqn.repl(tidx_caller);
            if isinf(njobs)
                callers_of_tidx_caller = find(lqn.taskgraph(:,tidx_caller));
                njobs = sum(mult(callers_of_tidx_caller)); %#ok<FNDSB>
                if isinf(njobs)
                    % if also the callers of tidx_caller are inf servers, then use
                    % an heuristic
                    njobs = min(sum(mult(isfinite(mult)) .* lqn.repl(isfinite(mult))),1e6);
                end
            end
            self.njobs(tidx_caller,idx) = njobs;
        else
            njobs = self.njobs(tidx_caller,idx);
        end
        caller_name = lqn.hashnames{tidx_caller};
        aidxClass{tidx_caller} = ClosedClass(model, caller_name, njobs, clientDelay);
        clientDelay.setService(aidxClass{tidx_caller}, Disabled.getInstance());
        for m=1:nreplicas
            serverStation{m}.setService(aidxClass{tidx_caller}, Disabled.getInstance());
        end
        aidxClass{tidx_caller}.completes = false;
        aidxClass{tidx_caller}.setReferenceClass(true); % renormalize residence times using the visits to the task
        aidxClass{tidx_caller}.attribute = [LayeredNetworkElement.TASK, tidx_caller];
        model.attribute.tasks(end+1,:) = [aidxClass{tidx_caller}.index, tidx_caller];
        %self.thinkproc
        clientDelay.setService(aidxClass{tidx_caller}, self.thinkproc{tidx_caller});
        if ~lqn.isref(tidx_caller)
            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, tidx_caller, 1, aidxClass{tidx_caller}.index];
        end
        for eidx = lqn.entriesof{tidx_caller}
            % create a class
            aidxClass{eidx} = ClosedClass(model, lqn.hashnames{eidx}, 0, clientDelay);
            clientDelay.setService(aidxClass{eidx}, Disabled.getInstance());
            for m=1:nreplicas
                serverStation{m}.setService(aidxClass{eidx}, Disabled.getInstance());
            end
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
            if isfield(lqn, 'arrival') && ~isempty(lqn.arrival) && ...
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
                    for m=1:nreplicas
                        serverStation{m}.setService(openClassForEntry, self.servtproc{bound_aidx});
                    end
                else
                    % Fallback to entry service (should not happen in well-formed models)
                    for m=1:nreplicas
                        serverStation{m}.setService(openClassForEntry, self.servtproc{eidx});
                    end
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
        if ishostlayer | any(any(lqn.issynccaller(tidx_caller, lqn.entriesof{idx}))) %#ok<OR2>
            % create a class
            aidxClass{aidx} = ClosedClass(model, lqn.hashnames{aidx}, 0, clientDelay);
            clientDelay.setService(aidxClass{aidx}, Disabled.getInstance());
            for m=1:nreplicas
                serverStation{m}.setService(aidxClass{aidx}, Disabled.getInstance());
            end
            aidxClass{aidx}.completes = false;
            aidxClass{aidx}.attribute = [LayeredNetworkElement.ACTIVITY, aidx];
            model.attribute.activities(end+1,:) = [aidxClass{aidx}.index, aidx];
            hidx = lqn.parent(lqn.parent(aidx)); % index of host processor
            if ~(ishostlayer && (hidx == idx))
                % set the host demand for the activity
                clientDelay.setService(aidxClass{aidx}, self.servtproc{aidx});
            end
            if lqn.sched(tidx_caller)~=SchedStrategy.REF % in 'ref' case the service activity is constant
                % updmap(end+1,:) = [idx, aidx, 1, idxClass{aidx}.index];
            end
            if iscachelayer && full(lqn.graph(eidx,aidx))
                clientDelay.setService(aidxClass{aidx}, self.servtproc{aidx});
            end

            % Create auxiliary think-time class if activity has think-time
            if ~isempty(lqn.actthink{aidx}) && lqn.actthink{aidx}.getMean() > GlobalConstants.FineTol
                aidxThinkClass{aidx} = ClosedClass(model, [lqn.hashnames{aidx},'.Think'], 0, clientDelay);
                aidxThinkClass{aidx}.completes = false;
                aidxThinkClass{aidx}.attribute = [LayeredNetworkElement.ACTIVITY, aidx];
                clientDelay.setService(aidxThinkClass{aidx}, lqn.actthink{aidx});
                for m=1:nreplicas
                    serverStation{m}.setService(aidxThinkClass{aidx}, Disabled.getInstance());
                end
                self.actthinkt_classes_updmap{idx}(end+1,:) = [idx, aidx, 1, aidxThinkClass{aidx}.index];
            end
        end
        % add a class for each outgoing call from this activity
        for cidx = lqn.callsof{aidx}
            callmean(cidx) = lqn.callproc{cidx}.getMean;
            switch lqn.calltype(cidx)
                case CallType.ASYNC
                    if lqn.parent(lqn.callpair(cidx,2)) == idx % add only if the target is serverStation
                        if ~hasSource % we need to add source and sink to the model
                            hasSource = true;
                            model.attribute.sourceIdx = length(model.nodes)+1;
                            sourceStation = Source(model,'Source');
                            sinkStation = Sink(model,'Sink');
                        end
                        cidxClass{cidx} = OpenClass(model, lqn.callhashnames{cidx}, 0);
                        sourceStation.setArrival(cidxClass{cidx}, Immediate.getInstance());
                        clientDelay.setService(cidxClass{cidx}, Disabled.getInstance());
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxClass{cidx}, Immediate.getInstance());
                        end
                        openClasses(end+1,:) = [cidxClass{cidx}.index, callmean(cidx), cidx];
                        model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                        cidxClass{cidx}.completes = false;
                        cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        minRespT = 0;
                        for tidx_act = lqn.actsof{idx}
                            minRespT = minRespT + lqn.hostdem{tidx_act}.getMean; % upper bound, uses all activities not just the ones reachable by this entry
                        end
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxClass{cidx}, Exp.fitMean(minRespT));
                        end
                    end
                case CallType.SYNC
                    cidxClass{cidx} = ClosedClass(model, lqn.callhashnames{cidx}, 0, clientDelay);
                    clientDelay.setService(cidxClass{cidx}, Disabled.getInstance());
                    for m=1:nreplicas
                        serverStation{m}.setService(cidxClass{cidx}, Disabled.getInstance());
                    end
                    model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                    cidxClass{cidx}.completes = false;
                    cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                    minRespT = 0;
                    for tidx_act = lqn.actsof{idx}
                        minRespT = minRespT + lqn.hostdem{tidx_act}.getMean; % upper bound, uses all activities not just the ones reachable by this entry
                    end
                    for m=1:nreplicas
                        serverStation{m}.setService(cidxClass{cidx}, Exp.fitMean(minRespT));
                    end
            end

            if callmean(cidx) ~= nreplicas
                switch lqn.calltype(cidx)
                    case CallType.SYNC
                        cidxAuxClass{cidx} = ClosedClass(model, [lqn.callhashnames{cidx},'.Aux'], 0, clientDelay);
                        cidxAuxClass{cidx}.completes = false;
                        cidxAuxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        clientDelay.setService(cidxAuxClass{cidx}, Immediate.getInstance());
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxAuxClass{cidx}, Disabled.getInstance());
                        end
                end
            end
        end
    end
end

% Ensure Source's sourceClasses and arrivalProcess arrays are properly sized for all classes
% This is needed because the Source may be created during class iteration
% when only some classes exist, and new closed classes added afterwards
% won't have corresponding entries in sourceClasses/arrivalProcess
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

P = model.initRoutingMatrix;
if hasSource
    for o = 1:size(openClasses,1)
        oidx = openClasses(o,1);
        p = 1 / openClasses(o,2); % divide by mean number of calls, they go to a server at random
        for m=1:nreplicas
            P{model.classes{oidx}, model.classes{oidx}}(sourceStation,serverStation{m}) = 1/nreplicas;
            for n=1:nreplicas
                P{model.classes{oidx}, model.classes{oidx}}(serverStation{m},serverStation{n}) = (1-p)/nreplicas;
            end
            P{model.classes{oidx}, model.classes{oidx}}(serverStation{m},sinkStation) = p;
        end
        cidx = openClasses(o,3); % 3 = source
        self.arvproc_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(sourceStation), oidx];
        for m=1:nreplicas
            self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), oidx];
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
    if ishostlayer
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
            end
        end
    end
    if (ishostlayer && hasDirectCallers) | any(any(lqn.issynccaller(tidx_caller, lqn.entriesof{idx}))) %#ok<OR2>
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
            P = recurActGraph(P, tidx_caller, eidx, aidxClass_eidx, jobPos);
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
        for m=1:nreplicas
            % Route: Source → ServerStation → Sink
            P{openClass, openClass}(sourceStation,serverStation{m}) = 1/nreplicas;
            P{openClass, openClass}(serverStation{m},sinkStation) = 1.0;
        end
    end
end

model.link(P);
self.ensemble{idx} = model;

    function [P, curClass, jobPos] = recurActGraph(P, tidx_caller, aidx, curClass, jobPos)
        jobPosKey(aidx) = jobPos;
        curClassKey{aidx} = curClass;
        nextaidxs = find(lqn.graph(aidx,:)); % these include the called entries
        if ~isempty(nextaidxs)
            isNextPrecFork(aidx) = any(isPostAndAct(nextaidxs)); % indexed on aidx to avoid losing it during the recursion
        end

        for nextaidx = nextaidxs % for all successor activities
            if ~isempty(nextaidx)
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
                            [P, jobPos, curClass] = routeSynchCall(P, jobPos, curClass);
                        case CallType.FWD
                            % Forwarding: after the source entry's activity completes,
                            % route to the forwarding target's bound activity on the same processor
                            if ishostlayer
                                target_entry = nextaidx;
                                target_bound_acts = find(lqn.graph(target_entry, :));
                                if ~isempty(target_bound_acts)
                                    target_aidx = target_bound_acts(1);  % First bound activity
                                    fwd_prob = full(lqn.graph(aidx, target_entry));  % Forwarding probability

                                    % Create class for target activity if needed
                                    if isempty(aidxClass{target_aidx})
                                        aidxClass{target_aidx} = ClosedClass(model, lqn.hashnames{target_aidx}, 0, clientDelay);
                                        clientDelay.setService(aidxClass{target_aidx}, Disabled.getInstance());
                                        for m=1:nreplicas
                                            serverStation{m}.setService(aidxClass{target_aidx}, Disabled.getInstance());
                                        end
                                        aidxClass{target_aidx}.completes = false;
                                        aidxClass{target_aidx}.attribute = [LayeredNetworkElement.ACTIVITY, target_aidx];
                                        model.attribute.activities(end+1,:) = [aidxClass{target_aidx}.index, target_aidx];
                                    end

                                    % Route from current position to target activity at server
                                    for m=1:nreplicas
                                        if jobPos == atClient
                                            P{curClass, aidxClass{target_aidx}}(clientDelay, serverStation{m}) = fwd_prob / nreplicas;
                                        else
                                            P{curClass, aidxClass{target_aidx}}(serverStation{m}, serverStation{m}) = fwd_prob;
                                        end
                                        serverStation{m}.setService(aidxClass{target_aidx}, lqn.hostdem{target_aidx});
                                    end
                                    self.servt_classes_updmap{idx}(end+1,:) = [idx, target_aidx, 2, aidxClass{target_aidx}.index];

                                    jobPos = atServer;
                                    curClass = aidxClass{target_aidx};

                                    % Check if target entry also has forwarding (chained forwarding)
                                    target_nextaidxs = find(lqn.graph(target_entry, :));
                                    for target_nextaidx = target_nextaidxs
                                        if target_nextaidx ~= target_aidx  % Skip the bound activity we just processed
                                            % Check if this is another forwarding call
                                            target_cidx = matchrow(lqn.callpair, [target_entry, target_nextaidx]);
                                            if ~isempty(target_cidx) && lqn.calltype(target_cidx) == CallType.FWD
                                                % Recursively handle chained forwarding
                                                [P, curClass, jobPos] = recurActGraph(P, tidx_caller, target_entry, curClass, jobPos);
                                            end
                                        end
                                    end
                                end
                            end
                    end
                else
                    % at this point, we have processed all calls, let us do the
                    % activities local to the task next
                    if isempty(intersect(lqn.eshift+(1:lqn.nentries), nextaidxs))
                        % if next activity is not an entry
                        jobPos = jobPosKey(aidx);
                        curClass = curClassKey{aidx};
                    else
                        if ismember(nextaidxs(find(nextaidxs==nextaidx)-1), lqn.eshift+(1:lqn.nentries))
                            curClassC = curClass;
                        end
                        jobPos = atClient;
                        curClass = curClassC;
                    end
                    if jobPos == atClient % at client node
                        if ishostlayer
                            if ~iscachelayer
                                for m=1:nreplicas
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(clientDelay, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, serverStation{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];
                                            P{curClass, forkClass}(clientDelay,joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,serverStation{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(clientDelay,serverStation{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    serverStation{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                    if isfunctionlayer
                                        serverStation{m}.setDelayOff(aidxClass{nextaidx}, lqn.setuptime{lqn.parent(nextaidx)}, lqn.delayofftime{lqn.parent(nextaidx)});
                                    end
                                end
                                jobPos = atServer;
                                curClass = aidxClass{nextaidx};
                                self.servt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 2, aidxClass{nextaidx}.index];
                            else
                                P{curClass, aidxClass{nextaidx}}(clientDelay,cacheNode) = full(lqn.graph(aidx,nextaidx));

                                cacheNode.setReadItemEntry(aidxClass{nextaidx},lqn.itemproc{aidx},lqn.nitems(aidx));
                                lqn.hitmissaidx = find(lqn.graph(nextaidx,:));
                                lqn.hitaidx = lqn.hitmissaidx(1);
                                lqn.missaidx = lqn.hitmissaidx(2);

                                cacheNode.setHitClass(aidxClass{nextaidx},aidxClass{lqn.hitaidx});
                                cacheNode.setMissClass(aidxClass{nextaidx},aidxClass{lqn.missaidx});

                                jobPos = atCache; % cache
                                curClass = aidxClass{nextaidx};
                                %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, lqn.hitaidx, 3, 3, aidxClass{nextaidx}.index, aidxClass{lqn.hitaidx}.index];
                                %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, lqn.missaidx, 3, 3, aidxClass{nextaidx}.index, aidxClass{lqn.missaidx}.index];
                            end
                        else % not ishostlayer
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
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.servtproc{nextaidx});
                            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    elseif jobPos == atServer || jobPos == atCache % at server station
                        if ishostlayer
                            if iscachelayer
                                curClass = aidxClass{nextaidx};
                                for m=1:nreplicas
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(cacheNode, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, serverStation{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];

                                            P{curClass, forkClass}(cacheNode,joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,serverStation{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(cacheNode,serverStation{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    serverStation{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                    %self.route_prob_updmap{idx}(end+1,:) = [idx, nextaidx, nextaidx, 3, 2, aidxClass{nextaidx}.index, aidxClass{nextaidx}.index];
                                end
                            else
                                for m=1:nreplicas
                                    if isNextPrecFork(aidx)
                                        % if next activity is a post-and
                                        P{curClass, curClass}(serverStation{m}, forkNode) = 1.0;
                                        f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                        forkClassStack(end+1) = curClass.index;
                                        P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                        P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, serverStation{m}) = 1.0;
                                    else
                                        if isPreAndAct(aidx)
                                            % before entering the job we go back to the entry class at the last fork
                                            forkClass = model.classes{forkClassStack(end)};
                                            forkClassStack(end) = [];
                                            P{curClass, forkClass}(serverStation{m},joinNode) = 1.0;
                                            P{forkClass, aidxClass{nextaidx}}(joinNode,serverStation{m}) = 1.0;
                                        else
                                            P{curClass, aidxClass{nextaidx}}(serverStation{m},serverStation{m}) = full(lqn.graph(aidx,nextaidx));
                                        end
                                    end
                                    serverStation{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                                end
                            end
                            jobPos = atServer;
                            curClass = aidxClass{nextaidx};
                            self.servt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 2, aidxClass{nextaidx}.index];
                        else
                            for m=1:nreplicas
                                if isNextPrecFork(aidx)
                                    % if next activity is a post-and
                                    P{curClass, curClass}(serverStation{m}, forkNode) = 1.0;
                                    f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                                    forkClassStack(end+1) = curClass.index;
                                    P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                                    P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, clientDelay) = 1.0;
                                else
                                    if isPreAndAct(aidx)
                                        % before entering the job we go back to the entry class at the last fork
                                        forkClass = model.classes{forkClassStack(end)};
                                        forkClassStack(end) = [];
                                        P{curClass, forkClass}(serverStation{m},joinNode) = 1.0;
                                        P{forkClass, aidxClass{nextaidx}}(joinNode,clientDelay) = 1.0;
                                    else
                                        P{curClass, aidxClass{nextaidx}}(serverStation{m},clientDelay) = full(lqn.graph(aidx,nextaidx));
                                    end
                                end
                            end
                            jobPos = atClient;
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.servtproc{nextaidx});
                            self.thinkt_classes_updmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    end
                    if aidx ~= nextaidx && ~isLoop
                        %% now recursively build the rest of the routing matrix graph
                        [P, curClass, jobPos] = recurActGraph(P, tidx_caller, nextaidx, curClass, jobPos);

                        % At this point curClassRec is the last class in the
                        % recursive branch, which we now close with a reply
                        if jobPos == atClient
                            P{curClass, aidxClass{tidx_caller}}(clientDelay,clientDelay) = 1;
                            if ~strcmp(curClass.name(end-3:end),'.Aux')
                                curClass.completes = true;
                            end
                        else
                            for m=1:nreplicas
                                P{curClass, aidxClass{tidx_caller}}(serverStation{m},clientDelay) = 1;
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

    function [P, jobPos, curClass] = routeSynchCall(P, jobPos, curClass)
        switch jobPos
            case atClient
                if lqn.parent(lqn.callpair(cidx,2)) == idx
                    % if a call to an entry of the server in this layer
                    if callmean(cidx) < nreplicas
                        P{curClass, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 - callmean(cidx); % note that callmean(cidx) < nreplicas
                        for m=1:nreplicas
                            %                                if isNextPrecFork(aidx)
                            %                                end
                            %                                     % if next activity is a post-and
                            %                                     P{curClass, curClass}(serverStation{m}, forkNode) = 1.0;
                            %                                     forkStackClass(end+1) = curClass.index;
                            %                                     f = find(nextaidx == nextaidxs(isPostAndAct(nextaidxs)));
                            %                                     P{curClass, curClass}(forkNode, forkOutputRouter{f}) = 1.0;
                            %                                     P{curClass, aidxClass{nextaidx}}(forkOutputRouter{f}, clientDelay) = 1.0;
                            %                                 else
                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = callmean(cidx) / nreplicas;
                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},clientDelay) = 1.0; % not needed, just to avoid leaving the Aux class disconnected
                        end
                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1.0; % not needed, just to avoid leaving the Aux class disconnected
                    elseif callmean(cidx) == nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1 / nreplicas;
                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},clientDelay) = 1.0;
                        end
                    else % callmean(cidx) > nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1 / nreplicas;
                            P{cidxClass{cidx}, cidxAuxClass{cidx}}(serverStation{m},clientDelay) = 1.0  ;
                            P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1.0 - 1.0 / (callmean(cidx) / nreplicas);
                        end
                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1.0 / (callmean(cidx));
                    end
                    jobPos = atClient;
                    clientDelay.setService(cidxClass{cidx}, Immediate.getInstance());
                    for m=1:nreplicas
                        serverStation{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                    end
                    curClass = cidxClass{cidx};
                else
                    % if it is not a call to an entry of the server
                    if callmean(cidx) < nreplicas
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = callmean(cidx)/nreplicas; % the mean number of calls is now embedded in the demand
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1; % the mean number of calls is now embedded in the demand
                        P{curClass, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 - callmean(cidx)/nreplicas; % the mean number of calls is now embedded in the demand
                        curClass = cidxAuxClass{cidx};
                    elseif callmean(cidx) == nreplicas
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > nreplicas
                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1; % the mean number of calls is now embedded in the demand
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;% / (callmean(cidx)/nreplicas); % the mean number of calls is now embedded in the demand
                        curClass = cidxAuxClass{cidx};
                    end
                    jobPos = atClient;
                    clientDelay.setService(cidxClass{cidx}, callservtproc{cidx});
                    self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                end
            case atServer % job at server
                if lqn.parent(lqn.callpair(cidx,2)) == idx
                    % if it is a call to an entry of the server
                    if callmean(cidx) < nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1 - callmean(cidx);
                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = callmean(cidx);
                            serverStation{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        end
                        jobPos = atClient;
                        curClass = cidxAuxClass{cidx};
                    elseif callmean(cidx) == nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1;
                        end
                        jobPos = atServer;
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1;
                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1 - 1 / (callmean(cidx));
                            P{cidxClass{cidx}, cidxAuxClass{cidx}}(serverStation{m},clientDelay) = 1 / (callmean(cidx));
                        end
                        jobPos = atClient;
                        curClass = cidxAuxClass{cidx};
                    end
                    for m=1:nreplicas
                        serverStation{m}.setService(cidxClass{cidx}, callservtproc{cidx});
                        self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                    end
                else
                    % if it is not a call to an entry of the server
                    % callmean not needed since we switched
                    % to ResidT to model service time at client
                    if callmean(cidx) < nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                        end
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxAuxClass{cidx};
                    elseif callmean(cidx) == nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                        end
                        curClass = cidxClass{cidx};
                    else % callmean(cidx) > nreplicas
                        for m=1:nreplicas
                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                        end
                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                        curClass = cidxAuxClass{cidx};
                    end
                    jobPos = atClient;
                    clientDelay.setService(cidxClass{cidx}, callservtproc{cidx});
                    self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                end
        end

        % After synch call returns, route through activity think-time class if applicable
        callingAidx = lqn.callpair(cidx, 1); % source activity of the call
        if ~isempty(aidxThinkClass{callingAidx})
            % Route from current class to think-time class at client
            P{curClass, aidxThinkClass{callingAidx}}(clientDelay, clientDelay) = 1.0;
            curClass = aidxThinkClass{callingAidx};
            jobPos = atClient;
        end
    end
end