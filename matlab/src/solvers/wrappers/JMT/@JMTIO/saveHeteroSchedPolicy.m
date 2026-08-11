function [simDoc, section] = saveHeteroSchedPolicy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEHETEROSCHEDPOLICY(SIMDOC, SECTION, IND)
% Saves heterogeneous scheduling policy to JMT XML.
% Generates the schedulingPolicy parameter.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
np = sn.nodeparam{ind};

% Check if this station has heterogeneous servers
if ~isfield(np, 'nservertypes') || np.nservertypes == 0
    return;
end

% Get scheduling policy for this station
if ~isfield(np, 'heteroschedpolicy')
    return;
end

policy = np.heteroschedpolicy;

% Create schedulingPolicy parameter
policyNode = simDoc.createElement('parameter');
policyNode.setAttribute('classPath', 'java.lang.String');
policyNode.setAttribute('name', 'schedulingPolicy');

valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(HeteroSchedPolicy.toJMTText(policy)));
policyNode.appendChild(valueNode);

section.appendChild(policyNode);
end
