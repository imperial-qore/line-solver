function [simElem, simDoc] = saveFCRMetrics(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEFCRMETRICS(SIMELEM, SIMDOC)
%
% Saves performance metrics for Finite Capacity Regions (blocking regions).
% Requests QLen, RespT, ResidT, and Tput metrics for each FCR.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.model.regions) || length(self.model.regions) == 0
    return;  % simElem, simDoc already assigned from inputs
end

% Metric types to request for FCRs
% JMT uses these names for blocking region metrics
metricTypes = {'Number of Customers', 'Response Time', 'Residence Time', 'Throughput'};
metricCounter = 0;

for r = 1:length(self.model.regions)
    fcrName = ['FCRegion', num2str(r)];

    for m = 1:length(metricTypes)
        metricType = metricTypes{m};
        performanceNode = simDoc.createElement('measure');
        performanceNode.setAttribute('alpha', num2str(1 - self.simConfInt, 2));
        performanceNode.setAttribute('name', ['FCR_', fcrName, '_', strrep(metricType, ' ', ''), '_', num2str(metricCounter)]);
        performanceNode.setAttribute('nodeType', 'region');
        performanceNode.setAttribute('precision', num2str(self.simMaxRelErr, 2));
        performanceNode.setAttribute('referenceNode', fcrName);
        performanceNode.setAttribute('referenceUserClass', '');  % FCR metrics are not class-specific
        performanceNode.setAttribute('type', metricType);
        performanceNode.setAttribute('verbose', 'false');
        simElem.appendChild(performanceNode);
        metricCounter = metricCounter + 1;
    end
end
end
