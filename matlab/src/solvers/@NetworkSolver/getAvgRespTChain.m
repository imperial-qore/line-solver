function [RN] = getAvgRespTChain(self,R)
% [RN] = GETAVGRESPTCHAIN(SELF,R)
% Return average response time aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[RNclass] = getAvgRespT(self);

[~,~,~,alpha] = sn_get_demands_chain(sn);

% compute average chain metrics
RN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    RN(:,c) = sum(RNclass(:,inchain).*alpha(:,inchain),2);
end

end
