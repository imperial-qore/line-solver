function [simDoc, section] = saveServersPerType(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVESERVERSPERTYPE(SIMDOC, SECTION, IND)
% Saves the number of servers per server type to JMT XML.
% Generates the serversPerServerType parameter array.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
ist = sn.nodeToStation(ind);

% Check if this station has heterogeneous servers
if isempty(sn.nservertypes) || sn.nservertypes(ist) == 0
    return;
end

% Get servers per type for this station
if isempty(sn.serverspertype) || ist > length(sn.serverspertype)
    return;
end

serversPerType = sn.serverspertype{ist};
if isempty(serversPerType)
    return;
end

% Create serversPerServerType parameter array
serversPerTypeNode = simDoc.createElement('parameter');
serversPerTypeNode.setAttribute('classPath', 'java.lang.Integer');
serversPerTypeNode.setAttribute('name', 'serversPerServerType');
serversPerTypeNode.setAttribute('array', 'true');

for t = 1:length(serversPerType)
    subNode = simDoc.createElement('subParameter');
    subNode.setAttribute('classPath', 'java.lang.Integer');
    subNode.setAttribute('name', 'serverTypesNumOfServers');

    valueNode = simDoc.createElement('value');
    valueNode.appendChild(simDoc.createTextNode(int2str(serversPerType(t))));
    subNode.appendChild(valueNode);
    serversPerTypeNode.appendChild(subNode);
end

section.appendChild(serversPerTypeNode);
end
