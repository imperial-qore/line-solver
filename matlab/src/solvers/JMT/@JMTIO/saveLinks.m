function [simElem, simDoc] = saveLinks(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVELINKS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
sn = self.getStruct;

[I,J]=find(sn.connmatrix);
for k=1:length(I)
    i=I(k);
    j=J(k);
    connectionNode = simDoc.createElement('connection');
    connectionNode.setAttribute('source', sn.nodenames(i));
    connectionNode.setAttribute('target', sn.nodenames(j));
    simElem.appendChild(connectionNode);
end
end
