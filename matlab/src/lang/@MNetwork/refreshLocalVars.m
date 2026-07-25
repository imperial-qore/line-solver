function nvars = refreshLocalVars(self)
% NVARS = REFRESHLOCALVARS()

R = self.getNumberOfClasses;
% Columns 1..R modulation phases, R+1..2R routing vars, 2R+1 the shared node
% block (cache width / BAS marker / polling controller). Columns 2R+1+r are the
% synchronous-call (REPLY) blocked-server counters, appended so that every
% existing nvars reader keeps its indices; see State.replyBlockInfo. They stay
% zero unless the model declares a REPLY signal, so no other model changes
% state width.
nvars = zeros(self.getNumberOfNodes, 3*R+1);
nodeparam = cell(self.getNumberOfNodes, 1);
rtnodes = self.sn.rtnodes;
% Draft SPN code:
% isp = [];
% ist = [];
% nodeToPlace = zeros(1, self.getNumberOfNodes);
% nodeToTransition = zeros(1, self.getNumberOfNodes);
% for ind=1:self.getNumberOfNodes
%     node = self.getNodeByIndex(ind);
%     switch class(node)
%     case 'Place'
%         isp = [isp, ind];
%         nodeToPlace(ind) = length(isp);
%     case 'Transition'
%         ist = [ist, ind];
%         nodeToTransition(ind) = length(ist);
%     end
% end
sn = self.sn;

for ind=1:self.getNumberOfNodes
    node = self.getNodeByIndex(ind);
    switch class(node)
        case 'Cache'
            nodeparam{ind} = struct();
            nodeparam{ind}.nitems = 0;
            nodeparam{ind}.accost = node.accessProb;
            for r=1:self.getNumberOfClasses
                if length(node.popularity) >= r && isa(node.popularity{r}, 'Distribution') && ~node.popularity{r}.isDisabled
                    nodeparam{ind}.nitems = max(nodeparam{ind}.nitems,node.popularity{r}.support(2));
                end
            end
            % nvars cache width = contents + (retrieval system) per-item occupancy bitmap
            if node.retrievalSystemCapacity > 0
                retrievalBitmapWidth = nodeparam{ind}.nitems;
            else
                retrievalBitmapWidth = 0;
            end
            nvars(ind,2*R+1) = node.totalCacheCapacity + retrievalBitmapWidth;
            nodeparam{ind}.itemcap = node.itemLevelCap;
            nodeparam{ind}.totalCacheCapacity = node.totalCacheCapacity;
            nodeparam{ind}.retrievalSystemCapacity = node.retrievalSystemCapacity;
            nodeparam{ind}.pread = cell(1,self.getNumberOfClasses);
            for r=1:self.getNumberOfClasses
                if length(node.popularity) < r || ~isa(node.popularity{r}, 'Distribution') || node.popularity{r}.isDisabled
                    nodeparam{ind}.pread{r} = NaN;
                else
                    nodeparam{ind}.pread{r} = node.popularity{r}.evalPMF(1:nodeparam{ind}.nitems);
                end
            end
            nodeparam{ind}.replacestrat = node.replacestrategy;
            if isprop(node,'admissionProb') && ~isempty(node.admissionProb)
                nodeparam{ind}.qlru = node.admissionProb;
            else
                nodeparam{ind}.qlru = 1.0;
            end
            % CLIMB rewritten as FIFO on unit-capacity lists with chain promotion -- see _kb/09-ldes-and-cache.md
            if node.replacestrategy == ReplacementStrategy.CLIMB
                Cclimb = sum(nodeparam{ind}.itemcap);
                nodeparam{ind}.itemcap = ones(1, Cclimb);
                nodeparam{ind}.replacestrat = ReplacementStrategy.FIFO;
                chainAc = diag(ones(1,Cclimb),1); chainAc(Cclimb+1,Cclimb+1) = 1;
                Kclimb = self.getNumberOfClasses;
                acClimb = cell(Kclimb, nodeparam{ind}.nitems);
                for vClimb = 1:Kclimb
                    for kClimb = 1:nodeparam{ind}.nitems
                        acClimb{vClimb,kClimb} = chainAc;
                    end
                end
                nodeparam{ind}.accost = acClimb;
            end
            nodeparam{ind}.hitclass = zeros(1,self.getNumberOfClasses);
            nodeparam{ind}.hitclass(1:length(node.server.hitClass)) = round(node.server.hitClass);
            nodeparam{ind}.missclass = zeros(1,self.getNumberOfClasses);
            nodeparam{ind}.missclass(1:length(node.server.missClass)) = round(node.server.missClass);

            % retrieval-system class matrix: pad to current class count and -1 outside known entries
            K = self.getNumberOfClasses;
            if isprop(node.server, 'retrievalClasses') && ~isempty(node.server.retrievalClasses)
                rc = -ones(nodeparam{ind}.nitems, K);
                rc(1:size(node.server.retrievalClasses,1), 1:min(size(node.server.retrievalClasses,2),K)) = ...
                    node.server.retrievalClasses(1:size(node.server.retrievalClasses,1), 1:min(size(node.server.retrievalClasses,2),K));
                nodeparam{ind}.retrievalClasses = rc;
            else
                nodeparam{ind}.retrievalClasses = -ones(max(nodeparam{ind}.nitems,1), K);
            end
            nodeparam{ind}.retrievalClassIndices = node.retrievalClassIndices;
            nodeparam{ind}.retrievalSystemQueueIndices = containers.Map('KeyType','int32','ValueType','any');
            keysList = keys(node.retrievalSystemQueueIndices);
            for kk = 1:numel(keysList)
                nodeparam{ind}.retrievalSystemQueueIndices(keysList{kk}) = node.retrievalSystemQueueIndices(keysList{kk});
            end

            % Store actual hit/miss/latency from solver results (if available)
            if isprop(node.server, 'actualHitProb') && ~isempty(node.server.actualHitProb)
                nodeparam{ind}.actualhitprob = full(node.server.actualHitProb);
                nodeparam{ind}.actualmissprob = full(node.server.actualMissProb);
            end
            if isprop(node.server, 'actualDelayedHitProb') && ~isempty(node.server.actualDelayedHitProb)
                nodeparam{ind}.actualdelayedhitprob = full(node.server.actualDelayedHitProb);
            end
            if isprop(node.server, 'actualResidT') && ~isempty(node.server.actualResidT)
                nodeparam{ind}.actualresidt = full(node.server.actualResidT);
            end
        case 'Fork'
            nodeparam{ind}.fanOut = node.output.tasksPerLink;
        case 'Join'
            nodeparam{ind}.joinStrategy = node.input.joinStrategy;
            nodeparam{ind}.fanIn = cell(1,self.getNumberOfClasses);
            for r=1:self.getNumberOfClasses
                nodeparam{ind}.fanIn{r} = node.input.joinRequired{r};
            end
        case 'Logger'
            nodeparam{ind}.fileName = node.fileName;
            nodeparam{ind}.filePath = node.filePath;
            nodeparam{ind}.startTime = node.getStartTime;
            nodeparam{ind}.loggerName = node.getLoggerName;
            nodeparam{ind}.timestamp = node.getTimestamp;
            nodeparam{ind}.jobID = node.getJobID;
            nodeparam{ind}.jobClass = node.getJobClass;
            nodeparam{ind}.timeSameClass = node.getTimeSameClass;
            nodeparam{ind}.timeAnyClass = node.getTimeAnyClass;            
        case {'Source'}
            for r=1:self.getNumberOfClasses
                if ~iscell(node.arrivalProcess) || length(node.arrivalProcess) < r || isempty(node.arrivalProcess{r})
                    continue;
                end
                switch class(node.arrivalProcess{r})
                    case {'MAP','MMPP2'} % Markov-modulated: track restart phase (see State.afterEvent ismkvmodclass)
                        nvars(ind,r) = nvars(ind,r) + 1;
                    case {'Replayer', 'Trace'}
                        if isempty(nodeparam{ind})
                            nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        end
                        if isempty(nodeparam{ind}{r})
                            nodeparam{ind}{r} = struct();
                        end
                        nodeparam{ind}{r}.(node.arrivalProcess{r}.params{1}.paramName) = node.arrivalProcess{r}.params{1}.paramValue;
                end
            end
        case {'Queue','QueueingStation','Delay','DelayStation','Transition'}
            for r=1:self.getNumberOfClasses
                switch class(node.server.serviceProcess{r}{3})
                    case {'MAP','MMPP2'} % Markov-modulated: track restart phase (see State.afterEvent ismkvmodclass)
                        nvars(ind,r) = nvars(ind,r) + 1;
                    case {'Replayer', 'Trace'}
                        if isempty(nodeparam{ind})
                            nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        end
                        if isempty(nodeparam{ind}{r})
                            nodeparam{ind}{r} = struct();
                        end
                        nodeparam{ind}{r}.(node.server.serviceProcess{r}{3}.params{1}.paramName) = node.server.serviceProcess{r}{3}.params{1}.paramValue;
                end
                if ~isa(node,'Delay') && ~isa(node,'Transition') && ~isempty(node.setupTime) && ~isempty(node.setupTime{r})
                    if isempty(nodeparam{ind})
                        nodeparam{ind} = cell(1,self.getNumberOfClasses);
                    end
                    if isempty(nodeparam{ind}{r})
                        nodeparam{ind}{r} = struct();
                    end
                    nodeparam{ind}{r}.setupTime = node.setupTime{r}.getProcess();
                    nodeparam{ind}{r}.delayoffTime = node.delayoffTime{r}.getProcess();
                end
                if (isa(node,'Queue') || isa(node,'QueueingStation'))
                    if isempty(nodeparam{ind}) && (~isempty(node.pollingType) || ~isempty(node.switchoverTime))
                        nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        for s=1:self.getNumberOfClasses
                            nodeparam{ind}{s} = struct();
                        end
                    end
                    if ~isempty(node.pollingType) && ~isempty(node.pollingType{r})
                        nodeparam{ind}{r}.pollingType = node.pollingType{r};
                        nodeparam{ind}{r}.pollingPar = node.pollingPar;
                    end
                    if ~isempty(node.switchoverTime) && ~isempty(node.switchoverTime{r})
                        if min(size(node.switchoverTime))==1
                            nodeparam{ind}{r}.switchoverTime = node.switchoverTime{r}.getProcess();
                            nodeparam{ind}{r}.switchoverProcId = ProcessType.toId(ProcessType.fromText(class(node.switchoverTime{r})));
                        else
                            for t=1:length(node.switchoverTime)
                                % Empty from-to switchover cells (diagonal, unset pairs) treated as Immediate, not dereferenced
                                so_rt = node.switchoverTime{r,t};
                                if isempty(so_rt)
                                    so_rt = Immediate();
                                end
                                nodeparam{ind}{r}.switchoverTime{t} = so_rt.getProcess();
                                nodeparam{ind}{r}.switchoverProcId(t) = ProcessType.toId(ProcessType.fromText(class(so_rt)));
                            end
                        end
                    end
                end
            end
            % Pass-and-swap (PAS) queues carry a node-level class
            % compatibility/swap graph rather than per-class parameters.
            if (isa(node,'Queue') || isa(node,'QueueingStation')) && ...
                    (SchedStrategy.toId(node.schedStrategy) == SchedStrategy.PAS || ...
                     SchedStrategy.toId(node.schedStrategy) == SchedStrategy.OI)
                if isempty(nodeparam{ind})
                    nodeparam{ind} = struct();
                end
                if SchedStrategy.toId(node.schedStrategy) == SchedStrategy.OI
                    % Order-independent queue: swap graph is always zero (empty),
                    % so class order is preserved on completion (plain OI).
                    nodeparam{ind}.swapGraph = zeros(R,R);
                elseif isempty(node.swapGraph)
                    % PAS default: complete compatibility graph (any class may
                    % swap with any other; no self-loops) when none was set.
                    nodeparam{ind}.swapGraph = ones(R,R) - eye(R);
                else
                    nodeparam{ind}.swapGraph = node.swapGraph;
                end
                % Total service rate function mu(c) of the ordered state vector c.
                nodeparam{ind}.svcRateFun = node.svcRateFun;
            end
    end

    for r=1:R
        switch sn.routing(ind,r)
            case RoutingStrategy.SQ
                nodeparam{ind}{r}.d = node.output.outputStrategy{r}{3}{1};
            case RoutingStrategy.WRROBIN
                nvars(ind,R+r) = nvars(ind,R+r) + 1;
                % save indexes of outgoing links
                if isempty(nodeparam) || isempty(nodeparam{ind}) % reinstantiate if not a cache
                    nodeparam{ind}{r} = struct();
                end
                nodeparam{ind}{r}.weights = zeros(1,self.sn.nnodes);
                nodeparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                for c=1:size(node.output.outputStrategy{1, r}{3},2)
                    destination = node.output.outputStrategy{1, r}{3}{c}{1};
                    weight = node.output.outputStrategy{1, r}{3}{c}{2};
                    nodeparam{ind}{r}.weights(destination.index) = weight;
                end
                % WRR cycle: replicate each outlink by its weight; afterEventRouter walks it one position per DEP (unit weights = plain RR)
                cycle = [];
                for ol = nodeparam{ind}{r}.outlinks
                    w = nodeparam{ind}{r}.weights(ol);
                    if w > 0
                        cycle = [cycle, repmat(ol, 1, max(1, round(w)))]; %#ok<AGROW>
                    else
                        cycle = [cycle, ol]; %#ok<AGROW>
                    end
                end
                nodeparam{ind}{r}.weighted_outlinks = cycle;
            case RoutingStrategy.RROBIN
                nvars(ind,R+r) = nvars(ind,R+r) + 1;
                % save indexes of outgoing links
                if isempty(nodeparam) || isempty(nodeparam{ind}) % reinstantiate if not a cache
                    nodeparam{ind}{r} = struct();
                end
                nodeparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
        end
    end
end

% True BAS blocking: shared trailing nvars column, isbasblocking/isbasdestination
% semantics and the BUG-83 upstream-vs-destination declaration forms -- see _kb/04-networkstruct.md (isbasblocking, isbasdestination fields)
isbasblocking = zeros(self.getNumberOfNodes, 1);
isbasdestnode = false(self.getNumberOfNodes, R);
for ind=1:self.getNumberOfNodes
    node = self.getNodeByIndex(ind);
    if isa(node,'Station') && ~isa(node,'Source') && ~isa(node,'Cache')
        [declares, destmask] = declaresBlockedMarker(self, ind, R);
        isbasdestnode = isbasdestnode | destmask;
        if declares
            % Option B (BUG-83): polling + BAS-blocking share one column, so a polling+BAS station is rejected here, before state-space build
            if self.sn.sched(self.sn.nodeToStation(ind)) == SchedStrategy.POLLING
                line_error(mfilename, sprintf(['True BAS blocking is not supported at a polling station (%s): the ', ...
                    'polling controller and the BAS blocked marker share one local-state column. Use a non-polling ', ...
                    'scheduling strategy at the blocking station, or remove the BAS drop rule.'], self.sn.nodenames{ind}));
            end
            nvars(ind, 2*R+1) = 1;
            isbasblocking(ind) = 1;
        end
    end
end

% Server breakdown status shares the trailing local-variable column with BAS/polling; combinations rejected before state-space build
hasbreakdown = zeros(self.getNumberOfNodes, 1);
if ~isempty(self.sn) && isfield(self.sn,'hasbreakdown') && ~isempty(self.sn.hasbreakdown)
    for ind=1:self.getNumberOfNodes
        if numel(self.sn.hasbreakdown) < ind || self.sn.hasbreakdown(ind) ~= 1
            continue
        end
        if isbasblocking(ind) == 1
            line_error(mfilename, sprintf(['Station ''%s'' combines server breakdowns with true-BAS blocking: the ' ...
                'breakdown status and the BAS blocked marker share one local-state column. Remove the BAS drop rule ' ...
                'or the breakdown.'], self.sn.nodenames{ind}));
        end
        if self.sn.sched(self.sn.nodeToStation(ind)) == SchedStrategy.POLLING
            line_error(mfilename, sprintf(['Server breakdowns are not supported at a polling station (%s): the polling ' ...
                'controller and the breakdown status share one local-state column.'], self.sn.nodenames{ind}));
        end
        nvars(ind, 2*R+1) = 1;
        hasbreakdown(ind) = 1;
    end
end

% Synchronous call (REPLY signal): a job of a class that expects a reply leaves
% this station for the callee but KEEPS its server, which is released only when
% the matching REPLY signal class arrives back here. Reserve one counter column
% per (station, calling class) so the state can carry the held servers; the job
% itself is at the callee and therefore absent from this station's marginal.
%
% The holding station is identified structurally, since a CTMC has no job
% identity to key on as LDES does: it is a station that the REPLY class is
% routed INTO. In the canonical shape Client -> Server -> (switch to Reply) ->
% Client, that selects the client and NOT the server -- the server's departure
% is the one that CREATES the reply, and LDES likewise does not block it
% (Solver_ssj: !classSwitchedToReply). Stations that never receive the reply
% class carry no counter and are untouched.
replyblock = zeros(self.getNumberOfNodes, R);
if isfield(self.sn,'syncreply') && ~isempty(self.sn.syncreply) && any(self.sn.syncreply >= 0)
    rtnodes = self.sn.rtnodes;
    nnodes = self.getNumberOfNodes;
    for r=1:R
        s = self.sn.syncreply(r) + 1; % stored 0-based (JAR convention)
        if s < 1 || s > R
            continue
        end
        for ind=1:nnodes
            if ~self.sn.isstation(ind) || self.sn.nodetype(ind) == NodeType.Source
                continue
            end
            ist = self.sn.nodeToStation(ind);
            % An INF station has a server for every job, so holding one is
            % immaterial and needs no state.
            if self.sn.sched(ist) == SchedStrategy.INF
                continue
            end
            % Does the reply class ever arrive here?
            arrivesHere = false;
            for i=1:nnodes
                for q=1:R
                    if rtnodes((i-1)*R+q, (ind-1)*R+s) > 0
                        arrivesHere = true;
                        break
                    end
                end
                if arrivesHere
                    break
                end
            end
            if ~arrivesHere
                continue
            end
            if SchedStrategy.toId(self.sn.sched(ist)) ~= SchedStrategy.FCFS
                line_error(mfilename, sprintf(['Synchronous calls (REPLY signals) are supported only at FCFS ', ...
                    'stations, but %s uses %s. Holding a server across a call has no representation in the ', ...
                    'state of the other disciplines.'], self.sn.nodenames{ind}, SchedStrategy.toText(self.sn.sched(ist))));
            end
            nvars(ind, 2*R+1+r) = 1;
            replyblock(ind, r) = 1;
        end
    end
end

% Polling controller width depends on discipline/switchover immediacy; State.pollingInfo is the single definition of the block layout
if ~isempty(self.sn)
    self.sn.nvars = nvars;
    self.sn.nodeparam = nodeparam;
    for ind=1:self.getNumberOfNodes
        if self.sn.isstation(ind) && self.sn.sched(self.sn.nodeToStation(ind)) == SchedStrategy.POLLING
            pinfo = State.pollingInfo(self.sn, ind);
            nvars(ind, 2*R+1) = pinfo.width;
            % memoize: State.pollingInfo is read once per state per
            % synchronization while the generator is built
            nodeparam{ind}{1}.pollinfo = pinfo;
            self.sn.nodeparam = nodeparam;
        end
    end
    self.sn.nvars = nvars;
    self.sn.isbasblocking = isbasblocking;
    % Project the node-indexed destination mask onto the station index space
    % used by sn.droprule, which is what State.arrivalIsLost indexes with.
    isbasdestination = false(self.sn.nstations, R);
    for ind=1:self.getNumberOfNodes
        if ~self.sn.isstation(ind)
            continue
        end
        jst = self.sn.nodeToStation(ind);
        if isnan(jst) || jst < 1
            continue
        end
        isbasdestination(jst,:) = isbasdestination(jst,:) | isbasdestnode(ind,:);
    end
    self.sn.isbasdestination = isbasdestination;
    self.sn.replyblock = replyblock;
end
end

function [tf, destmask] = declaresBlockedMarker(self, ind, R)
% [TF, DESTMASK] = DECLARESBLOCKEDMARKER(SELF, IND, R)
%
% DESTMASK is an (nnodes x R) logical marking the (destination node, class)
% pairs whose refusals must block IND rather than be lost.
%
% True when station IND is the BLOCKING (upstream) side of a true-BAS relation:
% it has a directly reachable destination station with a finite capacity, and
% BAS is declared EITHER on this station (upstream form, cqn_bas_blocking.m) OR
% on that destination (destination form, the JMT/LDES convention used by
% test_des_bas_closed.m; canonical since BUG-83). Both forms resolve to the same
% blocking station -- the held job sits here and the become-blocked edge starts
% here -- so the marker always lives on the upstream station. Whether the caller
% keys enumeration on the station's own drop rule (which fails for the
% destination form) or on the dedicated sn.isbasblocking field is what BUG-83
% fixes. Mirrors Network.declaresBlockedMarker in the JAR.
tf = false;
destmask = false(self.getNumberOfNodes, R);
node = self.getNodeByIndex(ind);
if ~isa(node,'Station') || isa(node,'Source') || isa(node,'Cache')
    return
end
dests = downstreamStations(self, ind);
for r = 1:R
    hereBAS = length(node.dropRule) >= r && node.dropRule(r) == DropStrategy.BAS;
    for d = 1:length(dests)
        dnode = self.getNodeByIndex(dests(d));
        if ~isa(dnode,'Station') || isa(dnode,'Source')
            continue
        end
        dcap = dnode.cap;
        if isempty(dcap) || ~isfinite(dcap) || dcap <= 0
            continue % the destination must be able to fill
        end
        thereBAS = length(dnode.dropRule) >= r && dnode.dropRule(r) == DropStrategy.BAS;
        if hereBAS || thereBAS
            % BAS declared upstream or on the reachable full destination. Do not
            % return early: every such destination must be recorded, since a
            % refusal at any of them has to block IND instead of dropping.
            tf = true;
            destmask(dests(d), r) = true;
        end
    end
end
end

function out = downstreamStations(self, ind)
% OUT = DOWNSTREAMSTATIONS(SELF, IND)
%
% Node indices of the stations directly downstream of IND, walking through
% intermediate stateless nodes (Router, ClassSwitch, ...) but stopping at the first
% station on each path, since a blocked job is held for its immediate destination.
out = [];
conn = self.sn.connmatrix;
if isempty(conn)
    return
end
n = self.getNumberOfNodes;
seen = false(1,n);
queue = ind;
seen(ind) = true;
while ~isempty(queue)
    cur = queue(1); queue(1) = [];
    for j = 1:min(n, size(conn,2))
        if conn(cur,j) ~= 1 || seen(j)
            continue
        end
        seen(j) = true;
        if isa(self.getNodeByIndex(j),'Station')
            out(end+1) = j; %#ok<AGROW>
        else
            queue(end+1) = j; %#ok<AGROW>
        end
    end
end
end
