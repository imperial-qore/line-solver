function P = makeErgodic(self, targetNode)
% P = MAKEERGODIC(TARGETNODE)
%
% Generates a modified routing matrix that makes the network ergodic
% by routing jobs from absorbing stations back to a target node.
%
% Input:
% TARGETNODE: (optional) Node or node name to route absorbing stations to.
%             If not specified, uses the first Delay node or first station.
%
% Output:
% P: Modified routing matrix (cell array) that can be used with model.link(P)
%
% Note: This method does NOT automatically relink the network. You must call
% model.link(P) with the returned routing matrix to apply the changes.
%
% Example:
%   P = model.makeErgodic();
%   model.link(P);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[isErg, info] = self.isRoutingErgodic();

if isErg
    line_warning(mfilename, 'Routing is already ergodic. No changes needed.');
    P = self.getLinkedRoutingMatrix();
    return;
end

if isempty(info.absorbingStations)
    line_warning(mfilename, 'No absorbing stations found, but network is not ergodic. Manual intervention required.');
    P = self.getLinkedRoutingMatrix();
    return;
end

% Determine target node
nodeNames = self.getNodeNames();
stationIdxs = self.getStationIndexes();

if nargin < 2 || isempty(targetNode)
    % Find first Delay node as default target
    targetIdx = [];
    for idx = stationIdxs
        if isa(self.nodes{idx}, 'Delay')
            targetIdx = idx;
            break;
        end
    end
    % If no Delay, use first station that's not absorbing
    if isempty(targetIdx)
        for idx = stationIdxs
            if ~ismember(nodeNames{idx}, info.absorbingStations)
                targetIdx = idx;
                break;
            end
        end
    end
    % Last resort: use first station
    if isempty(targetIdx) && ~isempty(stationIdxs)
        targetIdx = stationIdxs(1);
    end
    if isempty(targetIdx)
        line_error(mfilename, 'Cannot find a suitable target node for routing.');
    end
    targetName = nodeNames{targetIdx};
else
    if isa(targetNode, 'Node')
        targetName = targetNode.getName();
        targetIdx = self.getNodeIndex(targetName);
    else
        targetName = targetNode;
        targetIdx = self.getNodeIndex(targetName);
    end
    if isempty(targetIdx) || targetIdx <= 0
        line_error(mfilename, sprintf('Target node "%s" not found in network.', targetName));
    end
end

% Get current routing matrix
P = self.sn.rtorig;
if isempty(P)
    P = self.initRoutingMatrix();
end

K = self.getNumberOfClasses();

% Modify routing for each absorbing station
for i = 1:length(info.absorbingStations)
    absStationName = info.absorbingStations{i};
    absIdx = self.getNodeIndex(absStationName);

    if absIdx == targetIdx
        line_warning(mfilename, sprintf('Cannot route %s to itself. Skipping.', absStationName));
        continue;
    end

    % For each class, redirect self-loops to the target
    for r = 1:K
        for s = 1:K
            if ~isempty(P{r,s})
                % Clear all routing from this absorbing station
                P{r,s}(absIdx, :) = 0;
            end
        end
        % Route to target with same class
        if ~isempty(P{r,r})
            P{r,r}(absIdx, targetIdx) = 1.0;
        end
    end

    line_printf('Routing from %s redirected to %s for all classes.\n', absStationName, targetName);
end

line_printf('\nTo apply changes, call: model.link(P)\n');

end
