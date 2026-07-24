function Tard = getAvgTardHandles(self)
% TARD = GETAVGTARDHANDLES()
% Get handles for mean tardiness metrics.
%
% Tard(i,r): mean tardiness of class r at node i

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.handles) || ~isfield(self.handles,'Tard')
    % Initialize if needed
    self.getAvgHandles();
end

if isfield(self.handles,'Tard')
    Tard = self.handles.Tard;
else
    Tard = [];
end
end
