function [simDoc, section] = saveServerTypeNames(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVESERVERTYPENAMES(SIMDOC, SECTION, IND)
% Saves heterogeneous server type names to JMT XML.
% Generates the serverNames parameter array.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
ist = sn.nodeToStation(ind);

% Check if this station has heterogeneous servers
if isempty(sn.nservertypes) || sn.nservertypes(ist) == 0
    return;
end

% Get server type names for this station
if isempty(sn.servertypenames) || ist > length(sn.servertypenames)
    return;
end

names = sn.servertypenames{ist};
if isempty(names)
    return;
end

% Create serverNames parameter array
serverNamesNode = simDoc.createElement('parameter');
serverNamesNode.setAttribute('classPath', 'java.lang.String');
serverNamesNode.setAttribute('name', 'serverNames');
serverNamesNode.setAttribute('array', 'true');

for t = 1:length(names)
    subNode = simDoc.createElement('subParameter');
    subNode.setAttribute('classPath', 'java.lang.String');
    subNode.setAttribute('name', 'serverTypesNames');

    valueNode = simDoc.createElement('value');
    valueNode.appendChild(simDoc.createTextNode(names{t}));
    subNode.appendChild(valueNode);
    serverNamesNode.appendChild(subNode);
end

section.appendChild(serverNamesNode);
end
