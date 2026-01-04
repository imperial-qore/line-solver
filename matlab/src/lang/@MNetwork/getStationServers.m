function S = getStationServers(self)
% S = GETSTATIONSERVERS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

for i=1:getNumberOfStations(self)
    S(i,1) = self.stations{i}.numberOfServers;
end
end
