function info = getReducibilityInfo(self)
% INFO = GETREDUCIBILITYINFO()
%
% Returns detailed information about the reducibility structure of the
% network's routing matrix.
%
% Output:
% INFO: struct with the following fields:
%   - isRoutingErgodic: true if the routing is ergodic (irreducible)
%   - isReducible: true if the routing is reducible
%   - absorbingStations: cell array of names of absorbing stations
%   - transientStations: cell array of names of transient stations
%   - numSCCs: number of strongly connected components
%   - suggestedFixes: cell array of suggested routing fixes
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[isErg, info] = self.isRoutingErgodic();
info.isRoutingErgodic = isErg;
info.suggestedFixes = {};

if ~isErg
    % Generate suggested fixes for absorbing stations
    nodeNames = self.getNodeNames();
    stationIdxs = self.getStationIndexes();

    % Find a suitable return target (first Delay node, or first station)
    returnTarget = '';
    for idx = stationIdxs
        node = self.nodes{idx};
        if isa(node, 'Delay')
            returnTarget = nodeNames{idx};
            break;
        end
    end
    if isempty(returnTarget) && ~isempty(stationIdxs)
        returnTarget = nodeNames{stationIdxs(1)};
    end

    % Generate fix suggestions for each absorbing station
    for i = 1:length(info.absorbingStations)
        absStation = info.absorbingStations{i};
        if ~isempty(returnTarget) && ~strcmp(absStation, returnTarget)
            info.suggestedFixes{end+1} = sprintf(...
                'Route jobs from %s back to %s (e.g., P{class}(%s, %s) = 1.0)', ...
                absStation, returnTarget, absStation, returnTarget);
        end
    end
end

end
