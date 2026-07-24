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
cdscalingpeak = getLimitedClassDependencePeak(self);
jdscaling = getLimitedJointDependence(self);
jdscalingpeak = getLimitedJointDependencePeak(self);

%% init minimal structure
sn = NetworkStruct(); % create in self to ensure propagation
sn.isfjaugmented = self.isFJAugmented; % FJ tag-augmented copy marker (see ModelAdapter.fjtag)
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
% Detect Signal/OpenSignal/ClosedSignal classes; signaltarget semantics -- see _kb/04-networkstruct.md
sn.signaltarget = -ones(sn.nclasses, 1);
for r=1:sn.nclasses
    if isa(self.classes{r},'Signal') || isa(self.classes{r},'OpenSignal') || isa(self.classes{r},'ClosedSignal')
        sn.issignal(r) = true;
        sn.signaltype{r} = self.classes{r}.signalType;
        if ismethod(self.classes{r}, 'getTargetJobClassIndex')
            ti = self.classes{r}.getTargetJobClassIndex();
            if ~isempty(ti) && ti >= 1
                sn.signaltarget(r) = ti;
            end
        end
    end
end
% Initialize syncreply - maps each class to its expected reply signal class index (-1 if none)
sn.syncreply = -ones(sn.nclasses, 1);
for r=1:sn.nclasses
    if ~isempty(self.classes{r}.replySignalClass)
        sn.syncreply(r) = self.classes{r}.replySignalClass.index - 1; % 0-based for JAR
    end
end
% Initialize classspawn - class injected at the same station on each completion (-1 if none)
sn.classspawn = -ones(sn.nclasses, 1);
for r=1:sn.nclasses
    if isprop(self.classes{r}, 'spawnClass') && ~isempty(self.classes{r}.spawnClass)
        sn.classspawn(r) = self.classes{r}.spawnClass.index - 1; % 0-based for JAR
    end
end
% Initialize signal removal configuration fields
sn.signalremdist = cell(sn.nclasses,1);
sn.signalrempolicy = zeros(sn.nclasses, 1);
sn.iscatastrophe = false(sn.nclasses, 1);
% iscatastrophe derived from signaltype for every signal kind, not isCatastrophe() gated to open kinds -- see _kb/04-networkstruct.md
for r=1:sn.nclasses
    if isa(self.classes{r}, 'Signal') || isa(self.classes{r}, 'OpenSignal') ...
            || isa(self.classes{r}, 'ClosedSignal')
        if self.classes{r}.isCatastrophe()
            sn.iscatastrophe(r) = true;
        end
        sn.signalremdist{r} = self.classes{r}.removalDistribution;
        sn.signalrempolicy(r) = self.classes{r}.removalPolicy;
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
sn.cdscalingpeak = cdscalingpeak;
sn.jdscaling = jdscaling;
sn.jdscalingpeak = jdscalingpeak;
sn.nodetype = nodetypes;
sn.nstations = sum(sn.isstation);
sn.isstateful = (nodetypes == NodeType.Source | nodetypes == NodeType.Delay | nodetypes == NodeType.Queue | nodetypes == NodeType.Cache | nodetypes == NodeType.Join | nodetypes == NodeType.Router | nodetypes == NodeType.Place | nodetypes == NodeType.Transition | (nodetypes == NodeType.Fork & self.forkStateful));
sn.isstatedep = false(sn.nnodes,3); % col 1: buffer, col 2: srv, col 3: routing
sn.isfunction = [];
for ind = 1:sn.nstations
    if isa(self.stations{ind},'Queue')
        sn.isfunction(ind) = ~isempty(self.stations{ind}.setupTime);
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
            case {RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.RL, RoutingStrategy.SQ}
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
% sn.markidx populated BEFORE refreshProcesses so refreshProcessRepresentations can see the marked groups
if ~isempty(self.sn)
    self.sn.markidx = -ones(self.sn.nstations, self.sn.nclasses);
    for ist = 1:self.sn.nstations
        node = self.stations{ist};
        if isa(node, 'Source') && ~isempty(node.markedClasses)
            for k = 1:numel(node.markedClasses)
                self.sn.markidx(ist, node.markedClasses(k)) = k;
            end
        end
    end
end

refreshProcesses(self);

% Export patience/impatience fields for abandonment-aware solvers (MAPMsG, JMT, etc.)
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

% Initialize varsparam for cache item state tracking (mirrors JAR Network.java:4745-4747)
sn.varsparam = -ones(sn.nnodes, 1);

% Export balking and retrial fields for sn-driven solvers
sn.balkingStrategy = zeros(sn.nstations, sn.nclasses);
sn.balkingThresholds = cell(sn.nstations, sn.nclasses);
sn.retrialType = zeros(sn.nstations, sn.nclasses);
sn.retrialMu = zeros(sn.nstations, sn.nclasses);
sn.retrialPhi = zeros(sn.nstations, sn.nclasses);
sn.retrialProc = cell(sn.nstations, sn.nclasses);
sn.retrialMaxAttempts = -ones(sn.nstations, sn.nclasses);
sn.retrialPolicy = RetrialPolicy.LINEAR * ones(sn.nstations, sn.nclasses);
sn.orbitMaxJobs = -ones(sn.nstations, sn.nclasses);
sn.orbitImpatience = cell(sn.nstations, sn.nclasses);
sn.batchRejectProb = zeros(sn.nstations, sn.nclasses);
% Breakdown/repair rates are per STATION (server property); degraded service while down is per (station,class)
sn.hasbreakdown = zeros(sn.nnodes, 1);
sn.breakdownMu = zeros(sn.nstations, 1);
sn.repairMu = zeros(sn.nstations, 1);
sn.breakdownProc = cell(sn.nstations, 1);
sn.repairProc = cell(sn.nstations, 1);
sn.downServiceRates = zeros(sn.nstations, sn.nclasses);
for ist = 1:sn.nstations
    node = self.stations{ist};
    if isa(node, 'Queue') && ~isempty(node.breakdownFailure) && ~isempty(node.breakdownRepair)
        sn.hasbreakdown(sn.stationToNode(ist)) = 1;
        sn.breakdownMu(ist) = 1 / node.breakdownFailure.getMean();
        sn.repairMu(ist) = 1 / node.breakdownRepair.getMean();
        sn.breakdownProc{ist} = node.breakdownFailure.getProcess();
        sn.repairProc{ist} = node.breakdownRepair.getProcess();
        for r = 1:sn.nclasses
            dsvc = [];
            if numel(node.breakdownDownService) == 1
                dsvc = node.breakdownDownService{1};
            elseif numel(node.breakdownDownService) >= r
                dsvc = node.breakdownDownService{r};
            end
            if ~isempty(dsvc) && isa(dsvc, 'Distribution') && ~isa(dsvc, 'Disabled')
                if ~isa(dsvc, 'Exp')
                    line_error(mfilename, sprintf(['Station ''%s'': the down-server service distribution must be ' ...
                        'exponential; a phase-type degraded service would need its own phase block in the joint chain.'], ...
                        node.getName()));
                end
                dmean = dsvc.getMean();
                if dmean > 0
                    sn.downServiceRates(ist, r) = 1 / dmean;
                end
            end
        end
    end
end
for ist = 1:sn.nstations
    node = self.stations{ist};
    if isa(node, 'Queue')
        for r = 1:sn.nclasses
            % Balking
            [bStrat, bThresh] = node.getBalking(self.classes{r});
            if ~isempty(bStrat)
                if isnumeric(bStrat)
                    sn.balkingStrategy(ist, r) = bStrat;
                elseif isa(bStrat, 'BalkingStrategy') && isprop(bStrat, 'id')
                    sn.balkingStrategy(ist, r) = bStrat.id;
                else
                    sn.balkingStrategy(ist, r) = double(bStrat);
                end
                sn.balkingThresholds{ist, r} = bThresh;
            end
            % retrialPolicy/orbitMaxJobs default to the classical unbounded per-customer orbit (plain setRetrial)
            if length(node.retrialPolicies) >= r && node.retrialPolicies(r) > 0
                sn.retrialPolicy(ist, r) = node.retrialPolicies(r);
            end
            if length(node.orbitMaxJobs) >= r && ~isempty(node.orbitMaxJobs(r))
                sn.orbitMaxJobs(ist, r) = node.orbitMaxJobs(r);
            end
            % Retrial
            [rDist, rMax] = node.getRetrial(self.classes{r});
            if ~isempty(rDist) && ~isa(rDist, 'Disabled')
                if isprop(rDist, 'type')
                    sn.retrialType(ist, r) = rDist.type;
                elseif isa(rDist, 'Exp')
                    sn.retrialType(ist, r) = ProcessType.EXP;
                elseif isa(rDist, 'Erlang')
                    sn.retrialType(ist, r) = ProcessType.ERLANG;
                elseif isa(rDist, 'HyperExp')
                    sn.retrialType(ist, r) = ProcessType.HYPEREXP;
                elseif isa(rDist, 'Det')
                    sn.retrialType(ist, r) = ProcessType.DET;
                elseif isa(rDist, 'Gamma')
                    sn.retrialType(ist, r) = ProcessType.GAMMA;
                elseif isa(rDist, 'Pareto')
                    sn.retrialType(ist, r) = ProcessType.PARETO;
                elseif isa(rDist, 'Weibull')
                    sn.retrialType(ist, r) = ProcessType.WEIBULL;
                elseif isa(rDist, 'Lognormal')
                    sn.retrialType(ist, r) = ProcessType.LOGNORMAL;
                elseif isa(rDist, 'Uniform')
                    sn.retrialType(ist, r) = ProcessType.UNIFORM;
                elseif isa(rDist, 'Coxian')
                    sn.retrialType(ist, r) = ProcessType.COXIAN;
                elseif isa(rDist, 'APH')
                    sn.retrialType(ist, r) = ProcessType.APH;
                elseif isa(rDist, 'PH')
                    sn.retrialType(ist, r) = ProcessType.PH;
                else
                    sn.retrialType(ist, r) = ProcessType.PH;
                end
                if ismethod(rDist, 'getMean')
                    meanVal = rDist.getMean();
                    if meanVal > 0
                        sn.retrialMu(ist, r) = 1 / meanVal;
                    else
                        sn.retrialMu(ist, r) = Inf;
                    end
                end
                if ismethod(rDist, 'getSCV')
                    sn.retrialPhi(ist, r) = rDist.getSCV();
                else
                    sn.retrialPhi(ist, r) = 1.0;
                end
                if ismethod(rDist, 'getProcess')
                    sn.retrialProc{ist, r} = rDist.getProcess();
                end
                sn.retrialMaxAttempts(ist, r) = rMax;
            end
            % Orbit impatience (abandonment from the retrial orbit); store {D0,D1}
            oDist = node.getOrbitImpatience(self.classes{r});
            if ~isempty(oDist) && ~isa(oDist, 'Disabled') && ismethod(oDist, 'getProcess')
                sn.orbitImpatience{ist, r} = oDist.getProcess();
            end
            % Batch rejection probability (retrial queues: probability that an
            % arriving batch is rejected in full when it does not fit)
            sn.batchRejectProb(ist, r) = node.getBatchRejectProbability(self.classes{r});
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
                  SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO, ...
                  SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO, ...
                  SchedStrategy.LCFS, SchedStrategy.LCFSPR, ...
                  SchedStrategy.FCFSPR];
    if ~any(ismember(sn.sched, prioScheds))
        line_warning(mfilename, 'Priority classes are specified but no priority-aware scheduling policy is used in the model. Priorities will be ignored.');
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
% Bump the version so that anything cached against the previous compilation
% (see MNetwork.structVersion) is treated as stale.
self.structVersion = self.structVersion + 1;

if any(sn.fj(:)) && ~self.isFJAugmented % if there are forks
    % Recompute visits post-mmt: sum auxiliary-class visits into originals (skipped on FJ tag-augmented copies) -- see _kb/04-networkstruct.md
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
            % Pair via fjclassmap explicitly; positional pairing against sn.inchain{origChain} breaks under mmt confinement -- see _kb/04-networkstruct.md
            a = fsn.inchain{new_chain}(jaux);
            if a > length(fjclassmap) || fjclassmap(a) <= 0
                continue; % not an auxiliary class
            end
            j = fjclassmap(a); % original class mirrored by aux class a
            self.sn.nodevisits{origChain}(:,j) = sn.nodeparam{origFork}.fanOut*(X(:,j) + Vaux(:,jaux));
        end
    end
end

sn = refreshRegions(self);
self.sn = sn;

% Heterogeneous server fields populated LAST so refreshLocalVars (which rebuilds nodeparam) cannot wipe them
for ist = 1:self.sn.nstations
    if isa(self.stations{ist}, 'Queue') && ~isempty(self.stations{ist}.serverTypes)
        nodeIdx = self.sn.stationToNode(ist);
        nTypes = length(self.stations{ist}.serverTypes);
        self.sn.nodeparam{nodeIdx}.nservertypes = nTypes;
        self.sn.nodeparam{nodeIdx}.servertypenames = cell(1, nTypes);
        self.sn.nodeparam{nodeIdx}.serverspertype = zeros(1, nTypes);
        self.sn.nodeparam{nodeIdx}.servercompat = zeros(nTypes, self.sn.nclasses);
        % heterorates(t,r) = service rate of class r on a type-t server
        % (1/mean; 0 where incompatible or unset - falls back to base rate).
        self.sn.nodeparam{nodeIdx}.heterorates = zeros(nTypes, self.sn.nclasses);

        hsd = self.stations{ist}.heteroServiceDistributions;
        for t = 1:nTypes
            st = self.stations{ist}.serverTypes{t};
            self.sn.nodeparam{nodeIdx}.servertypenames{t} = st.getName();
            self.sn.nodeparam{nodeIdx}.serverspertype(t) = st.numOfServers;

            % Build compatibility matrix and per-(type,class) rates
            classMap = [];
            if ~isempty(hsd) && isKey(hsd, st.getName())
                classMap = hsd(st.getName());
            end
            for r = 1:self.sn.nclasses
                if st.isCompatible(self.classes{r})
                    self.sn.nodeparam{nodeIdx}.servercompat(t, r) = 1;
                    if ~isempty(classMap) && isKey(classMap, self.classes{r}.getName())
                        dist = classMap(self.classes{r}.getName());
                        mval = dist.getMean();
                        if mval > 0
                            self.sn.nodeparam{nodeIdx}.heterorates(t, r) = 1/mval;
                        end
                    end
                end
            end
        end

        % Get heterogeneous scheduling policy
        if ~isempty(self.stations{ist}.heteroSchedPolicy)
            self.sn.nodeparam{nodeIdx}.heteroschedpolicy = self.stations{ist}.heteroSchedPolicy;
        end
    end
end
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
