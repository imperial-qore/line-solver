function alpha = getLimitedLoadDependence(self)
% alpha = GETLIMITEDLOADDEPENDENCE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);

% First pass: find maximum size needed
maxsize = 0;
for ist=1:M
    mu = self.stations{ist}.lldScaling;
    if ~isempty(mu)
        maxsize = max(maxsize, length(mu));
    end
end

% Initialize with correct size and fill with ones
alpha = ones(M, maxsize);

% Second pass: populate with actual values
for ist=1:M
    mu = self.stations{ist}.lldScaling;
    if ~isempty(mu)
        alpha(ist, 1:length(mu)) = mu;
        % Replace any zeros with ones
        alpha(ist, alpha(ist,:)==0) = 1;
    end
end
end
