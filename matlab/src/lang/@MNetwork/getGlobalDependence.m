function phi = getGlobalDependence(self)
% phi = GETGLOBALDEPENDENCE()
%
% Returns the network-level globally state-dependent rate scaling handle
% declared through setGlobalDependence, or [] if the model has none. The handle
% receives the full (nstations x nclasses) per-class population matrix; see
% setGlobalDependence for the return-shape contract. ISEMPTY of this value is
% the test for "the model has no global dependence at all".

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

phi = self.gdScaling;
end
