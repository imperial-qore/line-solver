function [simDoc, section] = saveHeteroSchedPolicy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEHETEROSCHEDPOLICY(SIMDOC, SECTION, IND)
% Saves heterogeneous scheduling policy to JMT XML.
% Generates the schedulingPolicy parameter.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
ist = sn.nodeToStation(ind);

% Check if this station has heterogeneous servers
if isempty(sn.nservertypes) || sn.nservertypes(ist) == 0
    return;
end

% Get scheduling policy for this station
if isempty(sn.heteroschedpolicy) || ist > length(sn.heteroschedpolicy)
    return;
end

policy = sn.heteroschedpolicy(ist);

% Create schedulingPolicy parameter
policyNode = simDoc.createElement('parameter');
policyNode.setAttribute('classPath', 'java.lang.String');
policyNode.setAttribute('name', 'schedulingPolicy');

valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(HeteroSchedPolicy.toText(policy)));
policyNode.appendChild(valueNode);

section.appendChild(policyNode);
end
