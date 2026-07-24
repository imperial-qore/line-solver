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
                        distName = self.nodes{i}.server.serviceProcess{r}{3}.name;
                        if ~strcmp(distName,'Immediate') && ~strcmp(distName,'Disabled')
                            self.setUsedLangFeature(distName);
                        end
                        if self.nodes{i}.numberOfServers > 1
                            %self.setUsedLangFeature('MultiServer')
                        end
                        self.setUsedLangFeature(SchedStrategy.toFeature(self.nodes{i}.schedStrategy));
                        self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                    end
                case 'Router'
                    self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                case 'Source'
                    distName = self.nodes{i}.input.sourceClasses{r}{3}.name;
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
                case 'Join'
                    self.setUsedLangFeature('Join');
                    self.setUsedLangFeature('Joiner');
                case 'Sink'
                    self.setUsedLangFeature('Sink');
                case 'Cache'
                    self.setUsedLangFeature('CacheClassSwitcher');
                    self.setUsedLangFeature('Cache');
                    self.setUsedLangFeature(ReplacementStrategy.toFeature(self.nodes{i}.replacestrategy));
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
% Setup/delay-off times (setDelayOff): registered so that solvers ignoring
% them (CTMC, SSA, MVA, NC, FLD) reject the model instead of silently solving
% it setup-free.
for i=1:getNumberOfNodes(self)
    if isa(self.nodes{i},'Queue') && ~isempty(self.nodes{i}.setupTime)
        self.setUsedLangFeature('SetupDelayOff');
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
used = self.usedFeatures;
end
