function [RN] = getAvgNodeRespTChain(self,R)
% [RN] = GETAVGNODERESPTCHAIN(SELF,R)
% Return average response time aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[~,~,RNclass,~,~,~] = getAvgNode(self);

[~,~,~,alpha,~,~,~] = sn_get_demands_chain(sn);
% compute average chain metrics
RN = zeros(sn.nnodes, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nnodes
        if ~isempty(RNclass)
            if sn.isstation(i)
                RN(i,c) = RNclass(sn.nodeToStation(i),inchain)*alpha(sn.nodeToStation(i),inchain)';
            end
        end
    end
end
end
