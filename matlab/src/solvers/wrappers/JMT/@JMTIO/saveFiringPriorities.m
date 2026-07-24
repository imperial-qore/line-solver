function [simDoc, section] = saveFiringPriorities(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEFIRINGPRIORITIES(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

firingPrioritiesNode = simDoc.createElement('parameter');
firingPrioritiesNode.setAttribute('classPath', 'java.lang.Integer');
firingPrioritiesNode.setAttribute('name', 'firingPriorities');
firingPrioritiesNode.setAttribute('array', 'true');

sn = self.getStruct;
numOfModes = sn.nodeparam{ind}.nmodes;
for m=1:(numOfModes)
    
    subFiringPriorityNode = simDoc.createElement('subParameter');
    subFiringPriorityNode.setAttribute('classPath', 'java.lang.Integer');
    subFiringPriorityNode.setAttribute('name', 'firingPriority');
    
    valueNode = simDoc.createElement('value');
    firingPrio = sn.nodeparam{ind}.firingprio(m);
    
    if isinf(firingPrio)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        valueNode.appendChild(simDoc.createTextNode(int2str(firingPrio)));
    end
    
    subFiringPriorityNode.appendChild(valueNode);
    firingPrioritiesNode.appendChild(subFiringPriorityNode);
end

section.appendChild(firingPrioritiesNode);
end
