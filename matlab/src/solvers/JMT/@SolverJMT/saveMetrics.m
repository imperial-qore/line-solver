function [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEMETRICS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if isempty(self.handles)
    self.getAvgHandles;
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
end


