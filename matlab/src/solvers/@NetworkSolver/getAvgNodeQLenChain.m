function [QN] = getAvgNodeQLenChain(self,Q)
% [QN] = GETAVGNODEQLENCHAIN(SELF,Q)

% Return average queue-lengths aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[QNclass, ~] = getAvgNode(self);

% compute average chain metrics
if ~isempty(QNclass)
    QN = cell2mat(cellfun(@(chain) sum(QNclass(:,chain), 2), sn.inchain, 'UniformOutput', false));
else
    QN = zeros(sn.nnodes, sn.nchains);
end
end
