function [WN] = getAvgResidTChain(self,W)
% [WN] = GETAVGRESIDTCHAIN(SELF,W)
% Return average residence time aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[WNclass] = getAvgResidT(self);

% compute average chain metrics
WN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    WN(:,c) = sum(WNclass(:,inchain),2);
end
end
