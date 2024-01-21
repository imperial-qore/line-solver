function used = getUsedLangFeatures(self)
% USED = GETUSEDLANGFEATURES()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

self.initUsedFeatures;
if ~isempty(self.getIndexClosedClasses)
    self.setUsedFeatures('ClosedClass');
end
if ~isempty(self.getIndexOpenClasses)
    self.setUsedFeatures('OpenClass');
end

% Get attributes
for i=1:getNumberOfNodes(self)
    for r=1:getNumberOfClasses(self)
        try % not all nodes have all classes
            switch class(self.nodes{i})
                case {'Queue','QueueingStation','DelayStation','Delay'}
                    if ~isempty(self.nodes{i}.server.serviceProcess{r})
                        self.setUsedFeatures(self.nodes{i}.server.serviceProcess{r}{3}.name);
                        if self.nodes{i}.numberOfServers > 1
                            %self.setUsedFeatures('MultiServer')
                        end
                        self.setUsedFeatures(SchedStrategy.toFeature(self.nodes{i}.schedStrategy));
                        self.setUsedFeatures(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                    end
                case 'Router'
                    self.setUsedFeatures(RoutingStrategy.toFeature(self.nodes{i}.output.outputStrategy{r}{2}));
                case 'Source'
                    self.setUsedFeatures(self.nodes{i}.input.sourceClasses{r}{3}.name);
                    self.setUsedFeatures('Source');
                case 'Cache'
                    self.setUsedFeatures('CacheServer');
                    self.setUsedFeatures('CacheClassSwitcher');
                    self.setUsedFeatures('Cache');
                case 'ClassSwitch'
                    self.setUsedFeatures('StatelessClassSwitcher');
                    self.setUsedFeatures('ClassSwitch');
                case 'Fork'
                    self.setUsedFeatures('Fork');
                    self.setUsedFeatures('Forker');
                case 'Join'
                    self.setUsedFeatures('Join');
                    self.setUsedFeatures('Joiner');
                case 'Sink'
                    self.setUsedFeatures('Sink');
                case 'Cache'
                    self.setUsedFeatures('CacheClassSwitcher');
                    self.setUsedFeatures('Cache');
                case 'Transition'
                    self.setUsedFeatures('Transition');
                    self.setUsedFeatures('Enabling');
                    self.setUsedFeatures('Timing');
                    self.setUsedFeatures('Firing');
                case 'Place'
                    self.setUsedFeatures('Storage');
                    self.setUsedFeatures('Linkage');
                    self.setUsedFeatures('Place');
            end
        end
    end
end
used = self.usedFeatures;
end
