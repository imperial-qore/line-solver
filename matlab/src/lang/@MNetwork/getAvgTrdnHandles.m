function Trdn = getAvgTrdnHandles(self)
% TRDN = GETAVGTRDNHANDLES()
% Get handles for mean tardiness metrics.
%
% Trdn(i,r): mean tardiness of class r at node i

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.handles) || ~isfield(self.handles,'Trdn')
    % Initialize if needed
    self.getAvgHandles();
end

if isfield(self.handles,'Trdn')
    Trdn = self.handles.Trdn;
else
    Trdn = [];
end
end
