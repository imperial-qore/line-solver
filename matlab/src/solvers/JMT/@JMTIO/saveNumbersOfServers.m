function [simDoc, section] = saveNumbersOfServers(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVENUMBERSOFSERVERS(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

numbersOfServersNode = simDoc.createElement('parameter');
numbersOfServersNode.setAttribute('classPath', 'java.lang.Integer');
numbersOfServersNode.setAttribute('name', 'numbersOfServers');
numbersOfServersNode.setAttribute('array', 'true');

sn = self.getStruct;
numOfModes = sn.nodeparam{ind}.nmodes;
for m=1:(numOfModes)
    
    subNumberOfServersNode = simDoc.createElement('subParameter');
    subNumberOfServersNode.setAttribute('classPath', 'java.lang.Integer');
    subNumberOfServersNode.setAttribute('name', 'numberOfServers');
    
    valueNode = simDoc.createElement('value');
    
    nmodeservers = sn.nodeparam{ind}.nmodeservers(m);
    if isinf(nmodeservers)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        valueNode.appendChild(simDoc.createTextNode(int2str(nmodeservers)));
    end
    
    subNumberOfServersNode.appendChild(valueNode);
    numbersOfServersNode.appendChild(subNumberOfServersNode);
end

section.appendChild(numbersOfServersNode);
end
