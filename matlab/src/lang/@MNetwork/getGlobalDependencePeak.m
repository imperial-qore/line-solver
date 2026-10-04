function peak = getGlobalDependencePeak(self)
% peak = GETGLOBALDEPENDENCEPEAK()
%
% Returns an (nstations x nclasses) matrix of the declared peak (maximum)
% globally state-dependent rate scaling, used to normalize utilization as
% Util = T*S/peak (the same T*S/c convention as ordinary multiserver stations).
% Returns [] if the model declares no global dependence.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.gdScaling)
    peak = [];
    return
end
M = getNumberOfStations(self);
K = getNumberOfClasses(self);
peak = self.gdScalingPeak;
if isempty(peak)
    line_error(mfilename, 'The model declares a global dependence with no peak rate; use setGlobalDependence(phi, peakRate).');
end
if isscalar(peak)
    peak = peak * ones(M,K);
elseif isequal(size(peak),[M,1])
    peak = repmat(peak(:),1,K);
elseif ~isequal(size(peak),[M,K])
    line_error(mfilename, sprintf('peakRate must be a scalar, an (%d x 1) column or an (%d x %d) matrix.', M, M, K));
end
end
