function SysTrdn = getAvgSysTrdnHandles(self)
% SYSTRDN = GETAVGSYSTRDNHANDLES()
% Get handles for mean system tardiness metrics.
%
% SysTrdn(1,r): mean system tardiness for class r

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.handles) || ~isfield(self.handles,'SysTrdn')
    % Initialize if needed
    self.getAvgHandles();
end

if isfield(self.handles,'SysTrdn')
    SysTrdn = self.handles.SysTrdn;
else
    SysTrdn = [];
end
end
