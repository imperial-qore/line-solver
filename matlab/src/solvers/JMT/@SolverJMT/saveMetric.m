function [simElem, simDoc]=saveMetric(self, simElem, simDoc, handles)
for i=1:size(handles,1)
    for r=1:size(handles,2)
        currentPerformanceIndex = handles{i,r};
        if currentPerformanceIndex.disabled == 0
            performanceNode = simDoc.createElement('measure');
            performanceNode.setAttribute('alpha', num2str(1 - self.simConfInt,2));
            performanceNode.setAttribute('name', strcat('Performance_', int2str(i)));
            performanceNode.setAttribute('nodeType', 'station');
            performanceNode.setAttribute('precision', num2str(self.simMaxRelErr,2));
            if isempty(currentPerformanceIndex.station)
                performanceNode.setAttribute('referenceNode', '');
            else
                performanceNode.setAttribute('referenceNode', currentPerformanceIndex.station.name);
            end
            performanceNode.setAttribute('referenceUserClass', currentPerformanceIndex.class.name);
            performanceNode.setAttribute('type', currentPerformanceIndex.type);
            performanceNode.setAttribute('verbose', 'false');
            simElem.appendChild(performanceNode);
        end
    end
end
end