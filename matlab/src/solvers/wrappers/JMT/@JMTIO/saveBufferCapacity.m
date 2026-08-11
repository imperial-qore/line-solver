function [simDoc, section] = saveBufferCapacity(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEBUFFERCAPACITY(SIMDOC, SECTION, NODEIDX)
%
% LINE uses Kendall notation where cap = K = total system capacity
% (queue + in-service jobs). JMT's "size" parameter also represents
% total capacity K.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
sizeNode = simDoc.createElement('parameter');
sizeNode.setAttribute('classPath', 'java.lang.Integer');
sizeNode.setAttribute('name', 'size');
valueNode = simDoc.createElement('value');
ist = sn.nodeToStation(ind);
if ~sn.isstation(ind) || isinf(sn.cap(ist))
    valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
else
    if sn.cap(ist) == sum(sn.njobs)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        nservers = sn.nservers(ist);
        if isinf(nservers)
            % Infinite servers (delay node) - no buffer needed
            valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
        else
            % Send LINE's total capacity K directly to JMT
            valueNode.appendChild(simDoc.createTextNode(int2str(sn.cap(ist))));
        end
    end
end

sizeNode.appendChild(valueNode);
section.appendChild(sizeNode);
end
