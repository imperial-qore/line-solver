function [RN] = getAvgRespTChain(self,R)
% [RN] = GETAVGRESPTCHAIN(SELF,R)
% Return average response time aggregated by chain
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
%if nargin == 1
%    [Q] = getAvgHandles(self);
%end
[RNclass] = getAvgRespT(self);

% compute chain visits
alpha = zeros(sn.nstations,sn.nclasses);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        for k=inchain % for all classes within the chain ( a class belongs to a single chain, the reference station must be identical
            %                        for all classes within a chain )
            alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k)/sum(sn.visits{c}(sn.stationToStateful(sn.refstat(k)),inchain));
        end
    end
end
alpha(~isfinite(alpha))=0;

% compute average chain metrics
RN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        if ~isempty(RNclass)
            RN(i,c) = RNclass(i,inchain)*alpha(i,inchain)';
        end
    end
end
end
