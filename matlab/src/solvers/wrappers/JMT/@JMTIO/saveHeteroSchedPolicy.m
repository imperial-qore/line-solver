function [simDoc, section] = saveHeteroSchedPolicy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEHETEROSCHEDPOLICY(SIMDOC, SECTION, IND)
% Saves heterogeneous scheduling policy to JMT XML.
% Generates the schedulingPolicy parameter.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pools = self.serverPools(ind);
if isempty(pools)
    return;
end

% Create schedulingPolicy parameter
policyNode = simDoc.createElement('parameter');
policyNode.setAttribute('classPath', 'java.lang.String');
policyNode.setAttribute('name', 'schedulingPolicy');

valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(HeteroSchedPolicy.toJMTText(pools.policy)));
policyNode.appendChild(valueNode);

section.appendChild(policyNode);
end
