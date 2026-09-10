function used = getUsedLangFeatures(self)
% USED = GETUSEDLANGFEATURES()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.initUsedFeatures;
if ~isempty(self.getIndexClosedClasses)
    self.setUsedLangFeature('ClosedClass');
end
if ~isempty(self.getIndexOpenClasses)
    self.setUsedLangFeature('OpenClass');
end

% G-network signal classes. Signal is an unresolved placeholder that becomes
% OpenSignal or ClosedSignal at refreshStruct time, so all three names are
% inspected here: getUsedLangFeatures may run before resolution.
for r=1:getNumberOfClasses(self)
    jobclass = self.classes{r};
    if isa(jobclass,'Signal') || isa(jobclass,'OpenSignal') || isa(jobclass,'ClosedSignal')
        if isa(jobclass,'ClosedSignal')
            self.setUsedLangFeature('ClosedSignal');
        elseif isa(jobclass,'OpenSignal')
            self.setUsedLangFeature('OpenSignal');
        else % unresolved placeholder: resolves by presence of a Source
            if isempty(self.getIndexSourceNode)
                self.setUsedLangFeature('ClosedSignal');
            else
                self.setUsedLangFeature('OpenSignal');
            end
        end
        switch jobclass.signalType
            case SignalType.NEGATIVE
                self.setUsedLangFeature('SignalType_NEGATIVE');
            case SignalType.REPLY
                self.setUsedLangFeature('SignalType_REPLY');
            case SignalType.CATASTROPHE
                self.setUsedLangFeature('SignalType_CATASTROPHE');
        end
        % see _kb/04-networkstruct.md (node/process construction notes) for rationale
        if ~isempty(jobclass.removalDistribution)
            self.setUsedLangFeature('SignalBatchRemoval');
        end
        if ~isempty(jobclass.removalPolicy) && jobclass.removalPolicy ~= RemovalPolicy.RANDOM
            self.setUsedLangFeature('SignalRemovalPolicy');
        end
    end
end

% Finite capacity regions: solvers that ignore FCRs must be able to reject
% them through the registry, not only through the imperative checks in their
% own runAnalyzer. Trigger matches python get_used_lang_features exactly: at
% least one region on the model.
if ~isempty(self.regions)
    self.setUsedLangFeature('Region');
end

% Get attributes
for i=1:getNumberOfNodes(self)
    for r=1:getNumberOfClasses(self)
        try % not all nodes have all classes
            switch class(self.nodes{i})
                case {'Queue','QueueingStation','DelayStation','Delay'}
                    if ~isempty(self.nodes{i}.server.serviceProcess{r})
                        % Immediate/Disabled are internal placeholders, not
                        % user-facing distributions (same exclusion as JAR)
                        % getFeatureName, not name: it marks the most SPECIFIC
                        % registry entry (Cox2 for a two-phase Coxian, Trace for
                        % a Trace), which name cannot do without also moving the
                        % ProcessType. SolverFeatureSet.supports resolves an
                        % undeclared specialization against its generalization.
                        distName = self.nodes{i}.server.serviceProcess{r}{3}.getFeatureName();
                        if ~strcmp(distName,'Immediate') && ~strcmp(distName,'Disabled')
                            self.setUsedLangFeature(distName);
                        end
                        % A finite-server station serving several jobs at
                        % once. A Delay carries Inf here and is NOT one: the
                        % single-server recursions (aql, mvac, rqna, explicit,
                        % the aba..scb bounds, RCAT, the M/G/1 closed forms)
                        % answered a c-server station as one server of the
                        % same rate, and only a structural predicate could say
                        % so. Marked once per model; setUsedLangFeature is
                        % idempotent.
                        nserv = self.nodes{i}.numberOfServers;
                        if ~isempty(nserv) && isfinite(nserv) && nserv > 1
                            self.setUsedLangFeature('MultiServer');
                        end
                        self.setUsedLangFeature(SchedStrategy.toFeature(self.nodes{i}.schedStrategy));
                        self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                    end
                case 'Router'
                    self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                case 'Source'
                    distName = self.nodes{i}.input.sourceClasses{r}{3}.getFeatureName();
                    if ~strcmp(distName,'Immediate') && ~strcmp(distName,'Disabled')
                        self.setUsedLangFeature(distName);
                    end
                    self.setUsedLangFeature('Source');
                case 'ClassSwitch'
                    self.setUsedLangFeature('StatelessClassSwitcher');
                    self.setUsedLangFeature('ClassSwitch');
                case 'Fork'
                    self.setUsedLangFeature('Fork');
                    self.setUsedLangFeature('Forker');
                    % variable forking levels: a name declared but never marked
                    % is a name no solver can ever refuse, so mark all three here.
                    % self.nodes{i}, not a bare `node`: that name was never
                    % assigned in this function, and the enclosing TRY swallowed
                    % the error, so none of the three was ever marked.
                    if ~isempty(self.nodes{i}.output.tasksPerLinkByDest)
                        self.setUsedLangFeature('ForkFanoutVector');
                    end
                    if ~isempty(self.nodes{i}.output.tasksPerLinkDist)
                        self.setUsedLangFeature('ForkFanoutRandom');
                    end
                    if ~isempty(self.nodes{i}.output.branchProb)
                        self.setUsedLangFeature('ForkBranchProbability');
                    end
                case 'Join'
                    self.setUsedLangFeature('Join');
                    self.setUsedLangFeature('Joiner');
                case 'Sink'
                    self.setUsedLangFeature('Sink');
                case 'Cache'
                    self.setUsedLangFeature('CacheClassSwitcher');
                    self.setUsedLangFeature('Cache');
                    self.setUsedLangFeature(ReplacementStrategy.toFeature(self.nodes{i}.replacestrategy));
                    if isprop(self.nodes{i},'costCap') && ~isempty(self.nodes{i}.costCap)
                        % per-list storage cost caps with per-item sizes
                        self.setUsedLangFeature('CacheItemSize');
                    end
                    if ~isempty(self.nodes{i}.retrievalClassIndices)
                        % delayed-hit retrieval system (setRetrievalSystem):
                        % per-item retrieval classes; JMT cannot load these
                        self.setUsedLangFeature('CacheRetrieval');
                    end
                case 'Transition'
                    self.setUsedLangFeature('Transition');
                    self.setUsedLangFeature('Enabling');
                    self.setUsedLangFeature('Timing');
                    self.setUsedLangFeature('Firing');
                    % Inhibitor arcs: flag only when a finite threshold is set,
                    % so plain SPNs are not gated out of solvers lacking it.
                    if isprop(self.nodes{i},'inhibitingConditions')
                        for m=1:numel(self.nodes{i}.inhibitingConditions)
                            if any(isfinite(self.nodes{i}.inhibitingConditions{m}(:)))
                                self.setUsedLangFeature('Inhibiting');
                                break;
                            end
                        end
                    end
                case 'Place'
                    self.setUsedLangFeature('Storage');
                    self.setUsedLangFeature('Linkage');
                    self.setUsedLangFeature('Place');
                    if isprop(self.nodes{i},'queueing') && ~isempty(self.nodes{i}.queueing) && self.nodes{i}.queueing
                        self.setUsedLangFeature('QueueingPlace');
                    end
            end
        end
    end
end
% A self-looping class is a closed class whose routing returns to its own
% reference station (sn.isslc). It was declared by eight solvers and never
% marked, so a solver that omitted the name (the bounds, RCAT, QNS) was never
% asked; it is now refused unless it declares the name. Read off the class
% objects, since the placeholder classes above are.
for r=1:getNumberOfClasses(self)
    if isa(self.classes{r}, 'SelfLoopingClass')
        self.setUsedLangFeature('SelfLoopingClass');
        break
    end
end
% Batch arrivals (Source.setArrivalBatch): a batch releases several jobs at one
% arrival epoch, which changes the arrival stream itself. Only the LDES engine
% reads sn.arrivalbatch; every other solver would serve the single-arrival
% stream of the same rate, so the name is marked and left to LDES to declare.
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    if isa(nodeObj,'Source') && ~isempty(nodeObj.arrivalBatch) && ...
            any(~cellfun(@isempty, nodeObj.arrivalBatch(:)))
        self.setUsedLangFeature('BatchArrival');
        break
    end
end
% A retrial orbit (Queue.setRetrial / setOrbit): the same per-class test
% refreshStruct makes for sn.retrialProc and the JMT writer makes for the
% retrial Queue constructor, a configured delay that is not the Disabled
% placeholder. A solver that reads no sn.retrial* field would answer the
% model with the refused jobs simply lost, so the orbit is gated by name.
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    if ~isa(nodeObj,'Queue') || isempty(nodeObj.retrialDelays)
        continue
    end
    marked = false;
    for r=1:numel(nodeObj.retrialDelays)
        rDist = nodeObj.retrialDelays{r};
        if ~isempty(rDist) && ~isa(rDist,'Disabled')
            self.setUsedLangFeature('Retrial');
            marked = true;
            break
        end
    end
    if marked
        break
    end
end
% A station or per-class buffer that can BIND: the one predicate
% NetworkSolver.checkBindingCapacity refuses on (node-level caps against the
% class populations, open classes always bind, Cache models exempt), asked
% here so the refusal has a registry name and model.help can show it for a
% solver method that does not declare it.
if self.findBindingCapacity()
    self.setUsedLangFeature('FiniteCapacity');
end
% Limited load-dependent scaling (setLoadDependence): registered so that a
% solver which never reads sn.lldscaling rejects the model instead of returning
% the alpha == 1 answer. It was declared by SolverMVA/NC/CTMC/SSA/FLD/LDES and
% never emitted here, so the gate could not fire and SolverMAM and the fluid
% methods outside the closing family silently dropped the scaling.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Station') && ~isempty(self.nodes{i}.lldScaling)
        self.setUsedLangFeature('LoadDependence');
        break
    end
end
% Class-dependent scaling (setLimitedClassDependence): registered so that
% solvers ignoring the handle (e.g. JMT) reject the model instead of
% silently solving it as if the scaling were absent.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Station') && ~isempty(self.nodes{i}.lcdScaling)
        self.setUsedLangFeature('ClassDependence');
        break
    end
end
% Joint-dependent scaling (setJointDependence): the NON-product-form eta_i(n)
% case. Registered so that solvers not plumbing the handle reject the model
% instead of silently solving it as if the scaling were absent.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Station') && ~isempty(self.nodes{i}.ljdScaling)
        self.setUsedLangFeature('JointDependence');
        break
    end
end
% Globally state-dependent scaling (setGlobalDependence): phi(n) over the FULL
% network state, the Whittle-network primitive. Only SolverCTMC plumbs it, so
% every other solver must reject the model rather than solve it unscaled.
if ~isempty(self.gdScaling)
    self.setUsedLangFeature('GlobalDependence');
end
% Setup/delay-off times (setDelayOff): registered so that solvers ignoring
% them (CTMC, SSA, MVA, NC, FLD) reject the model instead of silently solving
% it setup-free.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Queue') && ~isempty(self.nodes{i}.setupTime)
        self.setUsedLangFeature('SetupDelayOff');
        break
    end
end
% Server parallelism (setServerParallelism): a job seizing n>1 servers changes
% the effective capacity of the station, so a solver that cannot honour it must
% reject the model rather than solve it as if every job seized one server.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Queue') && self.nodes{i}.hasServerParallelism()
        self.setUsedLangFeature('ServerParallelism');
        break
    end
end
% see _kb/04-networkstruct.md (node/process construction notes) for rationale
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    if isa(nodeObj,'Queue') && ~isempty(nodeObj.impatienceTypes)
        for r=1:numel(nodeObj.impatienceTypes)
            if ~isempty(nodeObj.impatienceTypes{r}) && nodeObj.impatienceTypes{r} == ImpatienceType.RENEGING
                self.setUsedLangFeature('Reneging');
                break
            end
        end
    end
end
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    if isa(nodeObj,'Queue') && ~isempty(nodeObj.balkingStrategies)
        for r=1:numel(nodeObj.balkingStrategies)
            if ~isempty(nodeObj.balkingStrategies{r})
                self.setUsedLangFeature('Balking');
                break
            end
        end
    end
end
% see _kb/04-networkstruct.md (node/process construction notes) for rationale
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    if isa(nodeObj,'Queue') && ~isempty(nodeObj.breakdownFailure) && ~isempty(nodeObj.breakdownRepair)
        self.setUsedLangFeature('Breakdown');
        break
    end
end
% A genuine QUORUM join (PARTIAL, k of n siblings with k < n) fires on the
% k-th branch completion, which is not their maximum, and the n-k stragglers
% are discarded on arrival: a solver that serves it as a full join returns a
% silent wrong answer, so the rule is gated by its own name. k >= n and
% k <= 0 are full joins and are not flagged.
if getNumberOfNodes(self) > 0
    conn = self.getConnectionMatrix;
    for i=1:getNumberOfNodes(self)
        nodeObj = self.nodes{i};
        if ~isa(nodeObj,'Join') || isempty(nodeObj.input) || ~isprop(nodeObj.input,'joinStrategy')
            continue
        end
        % siblings are counted at the FORK, as the engines do: its out-degree
        % times tasksPerLink. The join in-degree is the fallback.
        nsib = nnz(conn(:,i));
        joinof = nodeObj.joinOf;
        if ~isempty(joinof) && isobject(joinof) && joinof.index > 0 && joinof.index <= size(conn,1)
            w = 1;
            if isprop(joinof.output,'tasksPerLink') && ~isempty(joinof.output.tasksPerLink)
                w = max(1, round(joinof.output.tasksPerLink(1)));
            end
            nsib = nnz(conn(joinof.index,:)) * w;
        end
        for r=1:numel(nodeObj.input.joinStrategy)
            js = nodeObj.input.joinStrategy{r};
            if isempty(js) || js == JoinStrategy.STD
                continue
            end
            kreq = [];
            if r <= numel(nodeObj.input.joinRequired)
                kreq = nodeObj.input.joinRequired{r};
            end
            if ~isempty(kreq) && kreq > 0 && (nsib <= 0 || kreq < nsib)
                self.setUsedLangFeature('JoinPartial');
            end
        end
    end
end
% Heterogeneous server pools (Queue.addServerType): the station is served by
% several pools with their own counts, class compatibilities and per-(type,
% class) rates. A solver that reads only sn.nservers answers for a homogeneous
% station of the same total size, which is a different system, so the pools are
% gated by their own name rather than silently flattened.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Queue') && ~isempty(self.nodes{i}.serverTypes)
        self.setUsedLangFeature('HeteroServers');
        break
    end
end
% A FIFO depository (Place.setDepartureDiscipline) releases a served token only
% after the tokens that entered service before it, so which output transitions
% are enabled depends on the arrival order and not only on the marking. No
% solver implements it, hence a clean refusal instead of a NORMAL depository's
% answer under the user's name.
for i=1:getNumberOfNodes(self)
    nodeObj = self.nodes{i};
    % NORMAL is 0 and is also the value an unset slot carries (Place.m:135),
    % so a plain inequality covers both.
    if isa(nodeObj,'Place') && ~isempty(nodeObj.departureDiscipline) && ...
            any(nodeObj.departureDiscipline ~= DepartureDiscipline.NORMAL)
        self.setUsedLangFeature('DepartureDiscipline');
        break
    end
end
used = self.usedFeatures;
end
