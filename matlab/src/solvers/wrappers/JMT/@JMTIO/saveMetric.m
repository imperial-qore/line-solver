function [simElem, simDoc]=saveMetric(self, simElem, simDoc, handles)
sn = self.getStruct;

% Get exportable classes (handles cache classes and class-switching)
exportClasses = self.getExportableClasses();

for i=1:size(handles,1)
    for r=1:size(handles,2)
        currentPerformanceIndex = handles{i,r};
        if currentPerformanceIndex.disabled == 0
            % Skip classes that should not be exported to JMT
            classIdx = currentPerformanceIndex.class.index;
            if ~exportClasses(classIdx)
                continue;
            end

            performanceNode = simDoc.createElement('measure');
            performanceNode.setAttribute('alpha', num2str(1 - self.simConfInt,2));
            performanceNode.setAttribute('name', strcat('Performance_', int2str(i)));
            % System-level metrics (station is empty) use nodeType=""
            if isempty(currentPerformanceIndex.station)
                performanceNode.setAttribute('nodeType', '');
                performanceNode.setAttribute('referenceNode', '');
            else
                performanceNode.setAttribute('nodeType', 'station');
                performanceNode.setAttribute('referenceNode', currentPerformanceIndex.station.name);
            end
            performanceNode.setAttribute('precision', num2str(self.simMaxRelErr,2));
            performanceNode.setAttribute('referenceUserClass', currentPerformanceIndex.class.name);
            performanceNode.setAttribute('type', MetricType.toText(currentPerformanceIndex.type));
            performanceNode.setAttribute('verbose', 'false');
            simElem.appendChild(performanceNode);
        end
    end
end
end