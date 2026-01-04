function refreshStruct(self, hardRefresh)
% REFRESHSTRUCT()
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sanitize(self);
resolveSignals(self);  % Resolve Signal placeholders to OpenSignal or ClosedSignal
if nargin<2
    hardRefresh = true;
end

%% store invariant information
if self.hasStruct && ~hardRefresh
    rtorig = self.sn.rtorig; % this must be destroyed with resetNetwork
end

if self.hasStruct && ~hardRefresh
    nodetypes = sn.nodetypes;
    classnames = sn.classnames;
    nodenames = sn.nodenames;
    refstat = sn.refstat;
else
    nodetypes = getNodeTypes(self);
    classnames = getClassNames(self);
    nodenames = getNodeNames(self);
    refstat = getReferenceStations(self);

    % Append FCR names and types to node lists (only when refreshing)
    for f = 1:length(self.regions)
        fcr = self.regions{f};
        nodenames{end+1} = fcr.getName();
        nodetypes(end+1) = NodeType.Region;
    end
end
conn = self.getConnectionMatrix;
njobs = getNumberOfJobs(self);
numservers = getStationServers(self);
lldscaling = getLimitedLoadDependence(self);
cdscaling = getLimitedClassDependence(self);
[ljdscaling, ljdcutoffs] = getLimitedJointDependence(self);
[ljcdscaling, ljcdcutoffs] = getLimitedJointClassDependence(self);

%% init minimal structure
sn = NetworkStruct(); % create in self to ensure propagation
if isempty(self.sn)
    sn.rtorig = {};
    sn.reward = {};
else
    sn.rtorig = self.sn.rtorig;
    if isfield(self.sn, 'reward')
        sn.reward = self.sn.reward;  % preserve reward definitions
    else
        sn.reward = {};
    end
end
% sn.nnodes counts physical nodes only; FCRs are virtual nodes appended to nodenames/nodetypes
sn.nnodes = numel(self.nodes);
sn.nclasses = length(classnames);

%% get routing strategies
routing = zeros(sn.nnodes, sn.nclasses);
for ind=1:sn.nnodes
    for r=1:sn.nclasses
        if isempty(self.nodes{ind}.output.outputStrategy{r})
            routing(ind,r) = RoutingStrategy.DISABLED;
        else
            routing(ind,r) = RoutingStrategy.fromText(self.nodes{ind}.output.outputStrategy{r}{2});
        end
    end
end
sn.isslc = false(sn.nclasses,1);
for r=1:sn.nclasses
    if isa(self.classes{r},'SelfLoopingClass')
        sn.isslc(r) = true;
    end
end
sn.issignal = false(sn.nclasses,1);
sn.signaltype = cell(sn.nclasses,1);
for r=1:sn.nclasses
    sn.signaltype{r} = NaN;
end
% Detect Signal classes and populate issignal/signaltype
% Check for Signal, OpenSignal, and ClosedSignal classes
for r=1:sn.nclasses
    if isa(self.classes{r},'Signal') || isa(self.classes{r},'OpenSignal') || isa(self.classes{r},'ClosedSignal')
        sn.issignal(r) = true;
        sn.signaltype{r} = self.classes{r}.signalType;
    end
end
% Initialize syncreply - maps each class to its expected reply signal class index (-1 if none)
sn.syncreply = -ones(sn.nclasses, 1);
for r=1:sn.nclasses
    if ~isempty(self.classes{r}.replySignalClass)
        sn.syncreply(r) = self.classes{r}.replySignalClass.index - 1; % 0-based for JAR
    end
end
% Initialize signal removal configuration fields
sn.signalRemovalDist = cell(sn.nclasses,1);
sn.signalRemovalPolicy = zeros(sn.nclasses, 1);
sn.isCatastrophe = false(sn.nclasses, 1);
% Populate removal configuration from Signal and ClosedSignal classes
for r=1:sn.nclasses
    if isa(self.classes{r}, 'Signal') || isa(self.classes{r}, 'OpenSignal')
        if self.classes{r}.isCatastrophe()
            sn.isCatastrophe(r) = true;
        end
        sn.signalRemovalDist{r} = self.classes{r}.removalDistribution;
        sn.signalRemovalPolicy(r) = self.classes{r}.removalPolicy;
    elseif isa(self.classes{r}, 'ClosedSignal')
        sn.signalRemovalDist{r} = self.classes{r}.removalDistribution;
        sn.signalRemovalPolicy(r) = self.classes{r}.removalPolicy;
    end
end
sn.nclosedjobs = sum(njobs(isfinite(njobs)));
sn.nservers = numservers;
sn.isstation = (nodetypes == NodeType.Source | nodetypes == NodeType.Delay | nodetypes == NodeType.Queue | nodetypes == NodeType.Join | nodetypes == NodeType.Place);
sn.nstations = sum(sn.isstation);
sn.scv = ones(sn.nstations,sn.nclasses);
sn.njobs = njobs(:)';
sn.refstat = refstat;
sn.space = cell(sn.nstations,1);
sn.routing = routing;
sn.chains = [];
sn.lst = {};
sn.lldscaling = lldscaling;
sn.cdscaling = cdscaling;
sn.ljdscaling = ljdscaling;
sn.ljdcutoffs = ljdcutoffs;
sn.ljcdscaling = ljcdscaling;
sn.ljcdcutoffs = ljcdcutoffs;
sn.nodetype = nodetypes;
sn.nstations = sum(sn.isstation);
sn.isstateful = (nodetypes == NodeType.Source | nodetypes == NodeType.Delay | nodetypes == NodeType.Queue | nodetypes == NodeType.Cache | nodetypes == NodeType.Join | nodetypes == NodeType.Router | nodetypes == NodeType.Place | nodetypes == NodeType.Transition);
sn.isstatedep = false(sn.nnodes,3); % col 1: buffer, col 2: srv, col 3: routing
sn.isfunction = [];
for ind = 1:sn.nstations
    if isa(self.stations{ind},'Queue')
        sn.isfunction(ind) = ~isempty(self.stations{ind}.setupTime);
    end
end

% Populate heterogeneous server fields from Queue nodes
sn.nservertypes = zeros(sn.nstations, 1);
sn.servertypenames = cell(sn.nstations, 1);
sn.serverspertype = cell(sn.nstations, 1);
sn.servercompat = cell(sn.nstations, 1);
sn.heteroschedpolicy = zeros(sn.nstations, 1);
for ist = 1:sn.nstations
    if isa(self.stations{ist}, 'Queue') && ~isempty(self.stations{ist}.serverTypes)
        nTypes = length(self.stations{ist}.serverTypes);
        sn.nservertypes(ist) = nTypes;
        sn.servertypenames{ist} = cell(1, nTypes);
        sn.serverspertype{ist} = zeros(1, nTypes);
        sn.servercompat{ist} = zeros(nTypes, sn.nclasses);

        for t = 1:nTypes
            st = self.stations{ist}.serverTypes{t};
            sn.servertypenames{ist}{t} = st.getName();
            sn.serverspertype{ist}(t) = st.numOfServers;

            % Build compatibility matrix
            for r = 1:sn.nclasses
                if st.isCompatible(self.classes{r})
                    sn.servercompat{ist}(t, r) = 1;
                end
            end
        end

        % Get heterogeneous scheduling policy
        if ~isempty(self.stations{ist}.heteroSchedPolicy)
            sn.heteroschedpolicy(ist) = self.stations{ist}.heteroSchedPolicy;
        end
    end
end

for ind=1:sn.nnodes
    switch sn.nodetype(ind)
        case NodeType.Cache
            sn.isstatedep(ind,2) = true; % state dependent service
            %         case NodeType.Place
            %             self.nodes{ind}.init();
            %         case NodeType.Transition
            %             self.nodes{ind}.init(); % this erases enablingConditions
    end

    for r=1:sn.nclasses
        switch sn.routing(ind,r)
            case {RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.RL, RoutingStrategy.KCHOICES}
                sn.isstatedep(ind,3) = true; % state dependent routing
        end
    end
end

sn.nstateful = sum(sn.isstateful);
sn.state = cell(sn.nstations,1);
for i=1:sn.nstateful
    sn.state{i} = [];
end
sn.nodenames = nodenames;
sn.classnames = classnames;
sn.connmatrix = conn;

sn.nodeToStateful =[];
sn.nodeToStation =[];
sn.stationToNode =[];
sn.stationToStateful =[];
sn.statefulToNode =[];
sn.statefulToStation =[];
for ind=1:sn.nnodes
    sn.nodeToStateful(ind) = nd2sf(sn,ind);
    sn.nodeToStation(ind) = nd2st(sn,ind);
end
for ist=1:sn.nstations
    sn.stationToNode(ist) = st2nd(sn,ist);
    sn.stationToStateful(ist) = st2sf(sn,ist);
end
for isf=1:sn.nstateful
    sn.statefulToNode(isf) = sf2nd(sn,isf);
    sn.statefulToStation(isf) = sf2st(sn,isf);
end

% Populate immediate feedback matrix (station x class)
sn.immfeed = false(sn.nstations, sn.nclasses);
for ist=1:sn.nstations
    nodeIdx = sn.stationToNode(ist);
    node = self.nodes{nodeIdx};
    for r=1:sn.nclasses
        % Check station-level setting (Queue only)
        stationHas = false;
        if isa(node, 'Queue') && ~isempty(node.immediateFeedback)
            if ischar(node.immediateFeedback) && strcmp(node.immediateFeedback, 'all')
                stationHas = true;
            elseif iscell(node.immediateFeedback)
                stationHas = any(cellfun(@(x) x == r, node.immediateFeedback));
            end
        end
        % Check class-level setting
        classHas = self.classes{r}.immediateFeedback;
        sn.immfeed(ist, r) = stationHas || classHas;
    end
end

sn.fj = self.getForkJoins();
self.sn = sn;
refreshPriorities(self);
if exist('refreshDeadlines', 'file')
    refreshDeadlines(self);
else
    % Inline implementation if method not loaded
    K = getNumberOfClasses(self);
    classdeadline = zeros(1,K);
    for r=1:K
        classdeadline(r) = self.getClassByIndex(r).deadline;
    end
    if ~isempty(self.sn)
        self.sn.classdeadline = classdeadline;
    end
end
refreshProcesses(self);

% Export patience/impatience fields for abandonment-aware solvers (MAPMSG, JMT, etc.)
sn = self.sn;
sn.patienceProc = cell(sn.nstations, sn.nclasses);
sn.impatienceClass = zeros(sn.nstations, sn.nclasses);  % ImpatienceType (RENEGING, BALKING)
sn.impatienceType = zeros(sn.nstations, sn.nclasses);   % ProcessType (EXP, ERLANG, etc.)
sn.impatienceMu = zeros(sn.nstations, sn.nclasses);     % Rate parameter (1/mean)
sn.impatiencePhi = zeros(sn.nstations, sn.nclasses);    % SCV parameter
sn.impatiencePhases = zeros(sn.nstations, sn.nclasses);
sn.impatienceProc = cell(sn.nstations, sn.nclasses);
sn.impatiencePie = cell(sn.nstations, sn.nclasses);
for ist = 1:sn.nstations
    node = self.stations{ist};
    if isa(node, 'Queue')
        for r = 1:sn.nclasses
            patienceDist = node.getPatience(self.classes{r});
            if ~isempty(patienceDist) && ~isa(patienceDist, 'Disabled')
                % Convert patience distribution to MAP representation
                patienceMAP = patienceDist.getProcess();
                if iscell(patienceMAP) && length(patienceMAP) >= 2
                    sn.patienceProc{ist, r} = patienceMAP;
                    sn.impatienceProc{ist, r} = patienceMAP;

                    % Get the ProcessType of the patience distribution
                    if isprop(patienceDist, 'type')
                        sn.impatienceType(ist, r) = patienceDist.type;
                    elseif isa(patienceDist, 'Exp')
                        sn.impatienceType(ist, r) = ProcessType.EXP;
                    elseif isa(patienceDist, 'Erlang')
                        sn.impatienceType(ist, r) = ProcessType.ERLANG;
                    elseif isa(patienceDist, 'HyperExp')
                        sn.impatienceType(ist, r) = ProcessType.HYPEREXP;
                    elseif isa(patienceDist, 'Det')
                        sn.impatienceType(ist, r) = ProcessType.DET;
                    elseif isa(patienceDist, 'Gamma')
                        sn.impatienceType(ist, r) = ProcessType.GAMMA;
                    elseif isa(patienceDist, 'Pareto')
                        sn.impatienceType(ist, r) = ProcessType.PARETO;
                    elseif isa(patienceDist, 'Weibull')
                        sn.impatienceType(ist, r) = ProcessType.WEIBULL;
                    elseif isa(patienceDist, 'Lognormal')
                        sn.impatienceType(ist, r) = ProcessType.LOGNORMAL;
                    elseif isa(patienceDist, 'Uniform')
                        sn.impatienceType(ist, r) = ProcessType.UNIFORM;
                    elseif isa(patienceDist, 'Coxian')
                        sn.impatienceType(ist, r) = ProcessType.COXIAN;
                    elseif isa(patienceDist, 'APH')
                        sn.impatienceType(ist, r) = ProcessType.APH;
                    elseif isa(patienceDist, 'PH')
                        sn.impatienceType(ist, r) = ProcessType.PH;
                    else
                        % Default to PH for other Markovian distributions
                        sn.impatienceType(ist, r) = ProcessType.PH;
                    end

                    % Extract distribution parameters for JMT export
                    % Get rate (mu = 1/mean)
                    if ismethod(patienceDist, 'getMean')
                        meanVal = patienceDist.getMean();
                        if meanVal > 0
                            sn.impatienceMu(ist, r) = 1 / meanVal;
                        else
                            sn.impatienceMu(ist, r) = Inf;
                        end
                    else
                        % Estimate from MAP: mean = -pi * D0^(-1) * e
                        D0 = patienceMAP{1};
                        D1 = patienceMAP{2};
                        n = size(D0, 1);
                        pi = ones(1, n) / n;  % Approximate stationary distribution
                        meanVal = -pi * (D0 \ ones(n, 1));
                        sn.impatienceMu(ist, r) = 1 / meanVal;
                    end

                    % Get SCV (phi)
                    if ismethod(patienceDist, 'getSCV')
                        sn.impatiencePhi(ist, r) = patienceDist.getSCV();
                    else
                        sn.impatiencePhi(ist, r) = 1.0;  % Default to exponential SCV
                    end

                    % Get number of phases
                    if ismethod(patienceDist, 'getNumberOfPhases')
                        sn.impatiencePhases(ist, r) = patienceDist.getNumberOfPhases();
                    else
                        sn.impatiencePhases(ist, r) = size(patienceMAP{1}, 1);
                    end

                    % Get initial probability vector (pie)
                    n = size(patienceMAP{1}, 1);
                    D0 = patienceMAP{1};
                    D1 = patienceMAP{2};
                    % For PH: pie is from the MAP representation D1 = (-D0*e)*pie
                    exitRates = -D0 * ones(n, 1);
                    idx = find(exitRates > 1e-10, 1);
                    if ~isempty(idx)
                        sn.impatiencePie{ist, r} = D1(idx, :) / exitRates(idx);
                    else
                        sn.impatiencePie{ist, r} = ones(1, n) / n;
                    end
                end
                % Get impatience class (RENEGING, BALKING)
                impType = node.getImpatienceType(self.classes{r});
                if ~isempty(impType)
                    sn.impatienceClass(ist, r) = impType;
                end
            end
        end
    end
end
self.sn = sn;

% Check if priorities are specified but no priority-aware scheduling policy is used
sn = self.sn;
if ~all(sn.classprio == sn.classprio(1))
    % Priority classes exist, check if any station uses priority-aware scheduling
    prioScheds = [SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
                  SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO, ...
                  SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO];
    if ~any(ismember(sn.sched, prioScheds))
        line_warning(mfilename, 'Priority classes are specified but no priority-aware scheduling policy (PSPRIO, DPSPRIO, GPSPRIO, HOL, FCFSPRIO, LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO) is used in the model. Priorities will be ignored.');
    else
        % Display priority info unless silent
        global LINEVerbose;
        if isempty(LINEVerbose) || LINEVerbose ~= VerboseLevel.SILENT
            [minPrio, minIdx] = min(sn.classprio);
            [maxPrio, maxIdx] = max(sn.classprio);
            highestPrioClasses = find(sn.classprio == minPrio);
            lowestPrioClasses = find(sn.classprio == maxPrio);
            highNames = strjoin(arrayfun(@(i) sn.classnames{i}, highestPrioClasses, 'UniformOutput', false), ',');
            lowNames = strjoin(arrayfun(@(i) sn.classnames{i}, lowestPrioClasses, 'UniformOutput', false), ',');
            line_printf('Priority: highest=%s, lowest=%s\n', highNames, lowNames);
        end
    end
end

if any(nodetypes == NodeType.Cache)
    % this also refreshes the routing matrix and the visits
    refreshChains(self, false); % wantVisits
else
    % this also refreshes the routing matrix and the visits
    refreshChains(self, true); % wantVisits
end
sn = self.sn;
refclasses = getReferenceClasses(self);
refclass = zeros(1,sn.nchains);
for c=1:sn.nchains
    isect = intersect(sn.inchain{c},find(refclasses));
    if any(isect)
        refclass(c) = isect;
    end
end
sn.refclass = refclass;
self.sn = sn;
refreshLocalVars(self); % depends on chains (rtnodes)
refreshPetriNetNodes(self);
refreshSync(self); % this assumes that refreshChain is called before
refreshGlobalSync(self);

sn = self.sn;
self.hasStruct = true;

if any(sn.fj(:)) % if there are forks
    % try to recompute visits after mixed-model transformation:
    % we obtain the visits of auxiliary class from the transformed model
    % then we sum them to the original classes
    %try
    [nonfjmodel, fjclassmap, forkmap, fanOut] = ModelAdapter.mmt(self);
    if any(fanOut==1)
        % in this case, one of the forks is degenerate with a single
        % outgoing link so no visit correction is needed
        line_warning(mfilename,'The specified fork-join topology has partial support, only SolverJMT simulation results may be reliable.\n');
        %    return;
    end
    fsn = nonfjmodel.getStruct();
    for new_chain=(sn.nchains+1):fsn.nchains
        anyAuxClass = fsn.inchain{new_chain}(1);
        origFork = forkmap(anyAuxClass);
        origChain = find(sn.chains(:,fjclassmap(anyAuxClass))); % original chain of the class
        fsn.nodevisits{new_chain}(fsn.nodetype == NodeType.Source | fsn.nodetype == NodeType.Sink | fsn.nodetype == NodeType.Fork,:) = 0;
        Vaux = fsn.nodevisits{new_chain}(:,fsn.inchain{new_chain});
        if fsn.nnodes ~= sn.nnodes
            % Build mapping from fsn nodes to sn nodes by name
            % This handles cases where ClassSwitch/Source/Sink differ between models
            VauxMapped = zeros(sn.nnodes, size(Vaux,2));
            for fsnRow = 1:fsn.nnodes
                nodeName = fsn.nodenames{fsnRow};
                snRow = find(strcmp(sn.nodenames, nodeName), 1);
                if ~isempty(snRow)
                    % This fsn node exists in sn - copy its visit data
                    VauxMapped(snRow, :) = Vaux(fsnRow, :);
                end
                % Nodes not in sn (ClassSwitch, Source, Sink) are skipped
            end
            Vaux = VauxMapped;
        end
        X = sn.nodevisits{origChain};
        for jaux=1:length(fsn.inchain{new_chain})
            j = sn.inchain{origChain}(jaux); % class in the original chain
            self.sn.nodevisits{origChain}(:,j) = sn.nodeparam{origFork}.fanOut*(X(:,j) + Vaux(:,jaux));
        end
    end
end

sn = refreshRegions(self);
self.sn = sn;
end

function stat_idx = nd2st(sn, node_idx)
% STAT_IDX = ND2ST(NODE_IDX)

if sn.isstation(node_idx)
    stat_idx = at(cumsum(sn.isstation),node_idx);
else
    stat_idx = NaN;
end
end

function node_idx = st2nd(sn,stat_idx)
% NODE_IDX = ST2ND(SELF,STAT_IDX)

v = cumsum(sn.isstation) == stat_idx;
if any(v)
    node_idx =  find(v, 1);
else
    node_idx = NaN;
end
end

function sful_idx = st2sf(sn,stat_idx)
% SFUL_IDX = ST2SF(SELF,STAT_IDX)

sful_idx = nd2sf(sn,st2nd(sn,stat_idx));
end

function sful_idx = nd2sf(sn, node_idx)
% SFUL_IDX = ND2SF(NODE_IDX)

if sn.isstateful(node_idx)
    sful_idx = at(cumsum(sn.isstateful),node_idx);
else
    sful_idx = NaN;
end
end

function node_idx = sf2nd(sn,stat_idx)
% NODE_IDX = SF2ND(SELF,STAT_IDX)

v = cumsum(sn.isstateful) == stat_idx;
if any(v)
    node_idx =  find(v, 1);
else
    node_idx = NaN;
end
end

function stat_idx = sf2st(sn,sful_idx)
% STAT_IDX = SF2ST(SELF,SFUL_IDX)

stat_idx = nd2st(sn,sf2nd(sn,sful_idx));
end
