function peak = getLimitedClassDependencePeak(self)
% peak = GETLIMITEDCLASSDEPENDENCEPEAK()
%
% Returns an (nstations x nclasses) matrix of the declared peak (maximum)
% class-dependent rate scaling per class. Used to normalize utilization at
% class-dependent stations as Util = T*S/peak (the same T*S/c convention as
% ordinary multiserver stations). A station's scalar peak is broadcast across
% all classes; non class-dependent stations are left as NaN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
peak = NaN(M,K);
for i=1:M
    if ~isempty(self.stations{i}.lcdScaling)
        pk = self.stations{i}.lcdScalingPeak;
        if isempty(pk)
            line_error(mfilename, sprintf('Class-dependent station %d has no declared peak rate; use setClassDependence(beta, peakRatePerClass).', i));
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
