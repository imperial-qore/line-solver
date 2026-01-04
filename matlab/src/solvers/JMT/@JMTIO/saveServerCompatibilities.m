function [simDoc, section] = saveServerCompatibilities(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVESERVERCOMPATIBILITIES(SIMDOC, SECTION, IND)
% Saves server-class compatibility matrix to JMT XML.
% Generates the serverCompatibilities parameter array.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
ist = sn.nodeToStation(ind);

% Check if this station has heterogeneous servers
if isempty(sn.nservertypes) || sn.nservertypes(ist) == 0
    return;
end

% Get compatibility matrix for this station
if isempty(sn.servercompat) || ist > length(sn.servercompat)
    return;
end

compat = sn.servercompat{ist};
if isempty(compat)
    return;
end

[nTypes, nClasses] = size(compat);

% Create serverCompatibilities parameter array
compatNode = simDoc.createElement('parameter');
compatNode.setAttribute('classPath', 'java.lang.Object');
compatNode.setAttribute('name', 'serverCompatibilities');
compatNode.setAttribute('array', 'true');

% For each server type
for t = 1:nTypes
    typeNode = simDoc.createElement('subParameter');
    typeNode.setAttribute('classPath', 'java.lang.Boolean');
    typeNode.setAttribute('name', 'serverTypesCompatibilities');
    typeNode.setAttribute('array', 'true');

    % For each class
    for r = 1:nClasses
        classNode = simDoc.createElement('subParameter');
        classNode.setAttribute('classPath', 'java.lang.Boolean');
        classNode.setAttribute('name', 'compatibilities');

        valueNode = simDoc.createElement('value');
        if compat(t, r) > 0
            valueNode.appendChild(simDoc.createTextNode('true'));
        else
            valueNode.appendChild(simDoc.createTextNode('false'));
        end
        classNode.appendChild(valueNode);
        typeNode.appendChild(classNode);
    end
    compatNode.appendChild(typeNode);
end

section.appendChild(compatNode);
end
