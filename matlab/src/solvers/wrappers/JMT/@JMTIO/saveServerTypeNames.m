function [simDoc, section] = saveServerTypeNames(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVESERVERTYPENAMES(SIMDOC, SECTION, IND)
% Saves heterogeneous server type names to JMT XML.
% Generates the serverNames parameter array.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pools = self.serverPools(ind);
if isempty(pools)
    return;
end

names = pools.names;

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
