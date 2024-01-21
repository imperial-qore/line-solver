function alpha = getLimitedLoadDependence(self)
% alpha = GETLIMITEDLOADDEPENDENCE()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
alpha = zeros(M,0);
for i=1:M
    mu = self.stations{i}.lldScaling;
    alpha(i,1:length(mu)) = mu;    
end
alpha(alpha==0)=1;
end
