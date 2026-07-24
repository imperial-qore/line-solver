function peak = getLimitedJointDependencePeak(self)
% peak = GETLIMITEDJOINTDEPENDENCEPEAK()
%
% Returns an (nstations x nclasses) matrix of the declared peak (maximum)
% joint-dependent rate scaling per class. Used to normalize utilization at
% joint-dependent stations as Util = T*S/peak (the same T*S/c convention as
% ordinary multiserver stations). A station's scalar peak is broadcast across
% all classes; non joint-dependent stations are left as NaN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
peak = NaN(M,K);
for i=1:M
    if ~isempty(self.stations{i}.ljdScaling)
        pk = self.stations{i}.ljdScalingPeak;
        if isempty(pk)
            line_error(mfilename, sprintf('Class-dependent station %d has no declared peak rate; use setJointDependence(beta, peakRatePerClass).', i));
        end
        if isscalar(pk)
            peak(i,:) = pk;
        elseif numel(pk) == K
            peak(i,:) = pk(:)';
        else
            line_error(mfilename, sprintf('peakRatePerClass at station %d must be a scalar or a vector of length nclasses=%d.', i, K));
        end
    end
end
end
