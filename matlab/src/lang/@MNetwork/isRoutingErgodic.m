function [isErg, info] = isRoutingErgodic(self, P)
% [ISERG, INFO] = ISROUTINGERGODIC()
% [ISERG, INFO] = ISROUTINGERGODIC(P)
%
% Check if the queueing network routing matrix is ergodic (irreducible).
% This checks only the routing structure, not the full CTMC state space.
% A routing is ergodic if all stations communicate, meaning the routing
% matrix does not create absorbing states or disconnected components.
%
% Note: This is a necessary but not sufficient condition for the full
% Markov chain to be ergodic. State-dependent routing, finite buffers,
% or other factors may still cause reducibility in the actual CTMC.
%
% Input:
% P: (optional) routing matrix cell array. If not provided, it will be
%    retrieved via getLinkedRoutingMatrix(). Passing P directly avoids
%    calling getStruct() which can have side effects during link().
%
% Output:
% ISERG: true if the routing is ergodic, false otherwise
% INFO: struct with details about non-ergodicity:
%   - absorbingStations: cell array of station names that are absorbing
%   - transientStations: cell array of station names that are transient
%   - numSCCs: number of strongly connected components
%   - isReducible: true if the routing creates a reducible structure
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

info = struct();
info.absorbingStations = {};
info.transientStations = {};
info.numSCCs = 1;
info.isReducible = false;

% Get the routing matrix if not provided
if nargin < 2 || isempty(P)
    P = self.getLinkedRoutingMatrix();
end
if isempty(P)
    % No routing defined yet
    isErg = true;
    return;
end

K = self.getNumberOfClasses;
I = self.getNumberOfNodes;
nodeNames = self.getNodeNames;

% Build an aggregate adjacency matrix across all classes
% A transition exists from node i to node j if there's any class
% that can route from i to j
adjMatrix = zeros(I, I);
for r = 1:K
    for s = 1:K
        if ~isempty(P{r,s})
            % Handle routing matrices that may be smaller than I x I
            % (e.g., when user specifies routing only for stations)
            [nRows, nCols] = size(P{r,s});
            if nRows <= I && nCols <= I
                adjMatrix(1:nRows, 1:nCols) = adjMatrix(1:nRows, 1:nCols) + (P{r,s} > 0);
            else
                % Matrix is larger than expected (includes ClassSwitch nodes)
                adjMatrix = adjMatrix + (P{r,s}(1:I, 1:I) > 0);
            end
        end
    end
end
adjMatrix = adjMatrix > 0;  % Convert to logical

% Find strongly connected components
[scc, isrec] = stronglyconncomp(adjMatrix);
numSCCs = max(scc);
info.numSCCs = numSCCs;

% Check for absorbing stations (stations with only self-loops or no outgoing edges)
% An absorbing station is one where once entered, jobs never leave
stationIdxs = self.getStationIndexes();

for idx = stationIdxs
    % Check if this station has any outgoing transitions to other nodes
    outgoing = adjMatrix(idx, :);
    outgoing(idx) = 0;  % Exclude self-loop

    if sum(outgoing) == 0
        % No outgoing transitions except possibly self-loop
        % Check if there's routing defined for this node
        hasRouting = false;
        for r = 1:K
            for s = 1:K
                if ~isempty(P{r,s})
                    [nRows, ~] = size(P{r,s});
                    if idx <= nRows && any(P{r,s}(idx,:) > 0)
                        hasRouting = true;
                        break;
                    end
                end
            end
            if hasRouting
                break;
            end
        end

        if hasRouting
            % This station has routing but only to itself - it's absorbing
            info.absorbingStations{end+1} = nodeNames{idx};
        end
    end
end

% Identify transient stations (stations in transient SCCs)
for i = 1:numSCCs
    if ~isrec(i)
        % This SCC is transient
        sccNodes = find(scc == i);
        for nodeIdx = sccNodes
            if ismember(nodeIdx, stationIdxs)
                nodeName = nodeNames{nodeIdx};
                if ~ismember(nodeName, info.absorbingStations)
                    info.transientStations{end+1} = nodeName;
                end
            end
        end
    end
end

% The network is ergodic if:
% 1. There are no absorbing stations (other than Sink for open networks)
% 2. There is only one recurrent SCC that contains all stations
%    (excluding Source/Sink for open networks)

% Check for Sink node (which is legitimately absorbing in open networks)
sinkIdx = self.getIndexSinkNode();
sourceIdx = self.getIndexSourceNode();

% Remove Sink from absorbing list if present (it's expected to be absorbing)
if ~isempty(sinkIdx) && sinkIdx > 0
    sinkName = nodeNames{sinkIdx};
    info.absorbingStations = setdiff(info.absorbingStations, {sinkName});
end

% Count recurrent SCCs (excluding Source/Sink nodes)
recurrentSCCs = find(isrec);
numRecurrentStations = 0;
for i = recurrentSCCs
    sccNodes = find(scc == i);
    for nodeIdx = sccNodes
        if ismember(nodeIdx, stationIdxs)
            % Exclude source and sink
            if (isempty(sourceIdx) || nodeIdx ~= sourceIdx) && ...
               (isempty(sinkIdx) || nodeIdx ~= sinkIdx)
                numRecurrentStations = numRecurrentStations + 1;
            end
        end
    end
end

% Determine ergodicity
hasAbsorbingStations = ~isempty(info.absorbingStations);
hasMultipleRecurrentSCCs = sum(isrec) > 1;

% For closed networks: should have exactly 1 recurrent SCC containing all stations
% For open networks: Sink is allowed to be absorbing
if hasAbsorbingStations || hasMultipleRecurrentSCCs
    isErg = false;
    info.isReducible = true;
else
    isErg = true;
    info.isReducible = false;
end

end
