function [simDoc, section] = saveFiringWeights(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEFIRINGWEIGHTS(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

firingWeightsNode = simDoc.createElement('parameter');
firingWeightsNode.setAttribute('classPath', 'java.lang.Double');
firingWeightsNode.setAttribute('name', 'firingWeights');
firingWeightsNode.setAttribute('array', 'true');

sn = self.getStruct;
numOfModes = sn.nodeparam{ind}.nmodes;
for m=1:(numOfModes)
    
    subFiringWeightNode = simDoc.createElement('subParameter');
    subFiringWeightNode.setAttribute('classPath', 'java.lang.Double');
    subFiringWeightNode.setAttribute('name', 'firingWeight');
    
    valueNode = simDoc.createElement('value');
    firingWeights = sn.nodeparam{ind}.fireweight(m);
    
    if isinf(firingWeights)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        valueNode.appendChild(simDoc.createTextNode(int2str(firingWeights)));
    end
    
    subFiringWeightNode.appendChild(valueNode);
    firingWeightsNode.appendChild(subFiringWeightNode);
end

section.appendChild(firingWeightsNode);
end
