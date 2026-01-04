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

% Get attributes
for i=1:getNumberOfNodes(self)
    for r=1:getNumberOfClasses(self)
        try % not all nodes have all classes
            switch class(self.nodes{i})
                case {'Queue','QueueingStation','DelayStation','Delay'}
                    if ~isempty(self.nodes{i}.server.serviceProcess{r})
                        self.setUsedLangFeature(self.nodes{i}.server.serviceProcess{r}{3}.name);
                        if self.nodes{i}.numberOfServers > 1
                            %self.setUsedLangFeature('MultiServer')
                        end
                        self.setUsedLangFeature(SchedStrategy.toFeature(self.nodes{i}.schedStrategy));
                        self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                    end
                case 'Router'
                    self.setUsedLangFeature(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                case 'Source'
                    self.setUsedLangFeature(self.nodes{i}.input.sourceClasses{r}{3}.name);
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
                case 'Place'
                    self.setUsedLangFeature('Storage');
                    self.setUsedLangFeature('Linkage');
                    self.setUsedLangFeature('Place');
            end
        end
    end
end
used = self.usedFeatures;
end
