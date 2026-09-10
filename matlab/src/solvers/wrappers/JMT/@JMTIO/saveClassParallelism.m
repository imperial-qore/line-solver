function [simDoc, section] = saveClassParallelism(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVECLASSPARALLELISM(SIMDOC, SECTION, IND)
% Saves the per-class job parallelism to JMT XML.
% Generates the classParallelism parameter array (Server.serverNumRequired).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pools = self.serverPools(ind);
if isempty(pools)
    return;
end

sn = self.getStruct;

parNode = simDoc.createElement('parameter');
parNode.setAttribute('classPath', 'java.lang.Integer');
parNode.setAttribute('name', 'classParallelism');
parNode.setAttribute('array', 'true');

for r = 1:sn.nclasses
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    parNode.appendChild(refClassNode);

    subNode = simDoc.createElement('subParameter');
    subNode.setAttribute('classPath', 'java.lang.Integer');
    subNode.setAttribute('name', 'serverParallelism');

    valueNode = simDoc.createElement('value');
    valueNode.appendChild(simDoc.createTextNode(int2str(pools.parallelism(r))));
    subNode.appendChild(valueNode);
    parNode.appendChild(subNode);
end

section.appendChild(parNode);
end
