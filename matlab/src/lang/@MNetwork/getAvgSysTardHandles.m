function SysTard = getAvgSysTardHandles(self)
% SYSTARD = GETAVGSYSTARDHANDLES()
% Get handles for mean system tardiness metrics.
%
% SysTard(1,r): mean system tardiness for class r

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.handles) || ~isfield(self.handles,'SysTard')
    % Initialize if needed
    self.getAvgHandles();
end

if isfield(self.handles,'SysTard')
    SysTard = self.handles.SysTard;
else
    SysTard = [];
end
end
