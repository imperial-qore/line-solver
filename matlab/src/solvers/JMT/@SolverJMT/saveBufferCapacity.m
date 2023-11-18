function [simDoc, section] = saveBufferCapacity(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEBUFFERCAPACITY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sn = self.getStruct;
sizeNode = simDoc.createElement('parameter');
sizeNode.setAttribute('classPath', 'java.lang.Integer');
sizeNode.setAttribute('name', 'size');
valueNode = simDoc.createElement('value');
if ~sn.isstation(ind) || isinf(sn.cap(sn.nodeToStation(ind)))
    valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
else
    if sn.cap(sn.nodeToStation(ind)) == sum(sn.njobs)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        valueNode.appendChild(simDoc.createTextNode(int2str(sn.cap(sn.nodeToStation(ind)))));
    end
end

sizeNode.appendChild(valueNode);
section.appendChild(sizeNode);
end
