function rtTypes = getRoutingStrategies(self)
% RTTYPES = GETROUTINGSTRATEGIES()

rtTypes = zeros(self.getNumberOfNodes,self.getNumberOfClasses);
for ind=1:self.getNumberOfNodes
    for r=1:self.getNumberOfClasses
        RoutingStrategy.toId(self.nodes{ind}.output.outputStrategy{r}{2});
    end
end
end