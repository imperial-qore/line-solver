function cutoff = getGlobalDependenceCutoff(self)
% cutoff = GETGLOBALDEPENDENCECUTOFF()
%
% Returns the per-slot OPEN-class truncation declared through the third argument
% of setGlobalDependence, used only when the handle is materialized onto the JSON
% wire. Returns [] if the model declares no global dependence, and the default 10
% if one was declared without an explicit cutoff.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.gdScaling)
    cutoff = [];
    return
end
cutoff = self.gdScalingCutoff;
if isempty(cutoff)
    cutoff = 10;
end
end
