function [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEMETRICS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.handles)
    self.model.getAvgHandles();  % Create handles on model
    self.handles = self.model.handles;  % Copy to JMTIO
end
handles = self.handles;
[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.Q);
[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.U);
[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.R);
[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.T);
[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.A);
%[simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.W);
% JMT ResidT is inconsistently defined with LINE's on some
% difficult class switching cases, hence we recompute it at the
% level of the NetworkSolver class to preserve consistency.

% Save tardiness metrics
if isfield(handles,'Tard') && ~isempty(handles.Tard)
    [simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.Tard);
end
if isfield(handles,'SysTard') && ~isempty(handles.SysTard)
    [simElem, simDoc] = saveMetric(self, simElem, simDoc, handles.SysTard);
end

% Save FCR (blocking region) metrics
[simElem, simDoc] = saveFCRMetrics(self, simElem, simDoc);

% Save Cache Hit Rate metrics for Cache nodes
[simElem, simDoc] = saveCacheHitRateMetrics(self, simElem, simDoc);
end

function [simElem, simDoc] = saveCacheHitRateMetrics(self, simElem, simDoc)
% Request Cache Hit Rate metrics for each Cache node
sn = self.getStruct;
model = self.model;

for ind = 1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Cache
        nodeName = sn.nodenames{ind};
        hitclass = sn.nodeparam{ind}.hitclass;

        % Request Cache Hit Rate for each class that has a hit class defined
        for r = 1:length(hitclass)
            if hitclass(r) > 0
                hitClassName = sn.classnames{hitclass(r)};

                performanceNode = simDoc.createElement('measure');
                performanceNode.setAttribute('alpha', num2str(1 - self.simConfInt, 2));
                performanceNode.setAttribute('name', sprintf('CacheHitRate_%s_%s', nodeName, hitClassName));
                performanceNode.setAttribute('nodeType', 'station');
                performanceNode.setAttribute('precision', num2str(self.simMaxRelErr, 2));
                performanceNode.setAttribute('referenceNode', nodeName);
                performanceNode.setAttribute('referenceUserClass', hitClassName);
                performanceNode.setAttribute('type', 'Cache Hit Rate');
                performanceNode.setAttribute('verbose', 'false');
                simElem.appendChild(performanceNode);
            end
        end
    end
end
end


