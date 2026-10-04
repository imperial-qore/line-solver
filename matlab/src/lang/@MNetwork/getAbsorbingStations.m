function [absorbingStations, absorbingIdxs] = getAbsorbingStations(self)
% [ABSORBINGSTATIONS, ABSORBINGIDXS] = GETABSORBINGSTATIONS()
%
% Returns the stations that are absorbing (i.e., once a job enters,
% it never leaves). These create reducible Markov chains.
%
% Output:
% ABSORBINGSTATIONS: cell array of Station objects that are absorbing
% ABSORBINGIDXS: vector of node indices of absorbing stations
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[~, info] = self.isRoutingErgodic();

absorbingStations = {};
absorbingIdxs = [];

for i = 1:length(info.absorbingStations)
    stationName = info.absorbingStations{i};
    station = self.getStationByName(stationName);
    if isa(station, 'Node')
        absorbingStations{end+1} = station;
        absorbingIdxs(end+1) = self.getNodeIndex(stationName);
    end
end

end
