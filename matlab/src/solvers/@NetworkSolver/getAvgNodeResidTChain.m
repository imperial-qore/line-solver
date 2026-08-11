function [WN] = getAvgNodeResidTChain(self,R)
% [WN] = GETAVGNODERESIDTCHAIN(SELF,R)
% Return average response time aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[~,~,~,~,~,WNclass] = getAvgNode(self);

[~,~,~,alpha,~,~,~] = sn_get_demands_chain(sn);
% compute average chain metrics
WN = zeros(sn.nnodes, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nnodes
        if ~isempty(WNclass)
            if sn.isstation(i)
                WN(i,c) = WNclass(sn.nodeToStation(i),inchain)*alpha(sn.nodeToStation(i),inchain)';
            end
        end
    end
end
end
