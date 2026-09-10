function [simDoc, section] = saveNumberOfServers(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVENUMBEROFSERVERS(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
sizeNode = simDoc.createElement('parameter');
sizeNode.setAttribute('classPath', 'java.lang.Integer');
sizeNode.setAttribute('name', 'maxJobs');

sn = self.getStruct;
ist = sn.nodeToStation(ind);

% For LPS, use 1 server (single server with PS). The LPS limit is handled
% by creating an implicit FCR region in saveRegions.
if sn.sched(ist) == SchedStrategy.LPS
    maxJobs = 1;  % LPS uses single server with FCR for admission control
else
    % Use the maximum of nservers and lldscaling (for load-dependent models).
    % Load-dependent scaling with min(1:N,c) represents a c-server queue,
    % where max(lldscaling) = c.
    maxJobs = sn.nservers(ist);
    if ~isempty(sn.lldscaling) && size(sn.lldscaling, 1) >= ist
        maxJobs = max(maxJobs, max(sn.lldscaling(ist,:)));
    end
end

valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(maxJobs)));

sizeNode.appendChild(valueNode);
section.appendChild(sizeNode);
end
