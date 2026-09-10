function [AN] = getAvgNodeArvRChain(self,A)
% [AN] = GETAVGNODEARVRCHAIN(SELF,A)
% Return average arrival rates aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[~,~,~,~,ANclass, ~] = getAvgNode(self);

% compute average chain metrics
if ~isempty(ANclass)
    AN = cell2mat(cellfun(@(chain) sum(ANclass(:,chain), 2), sn.inchain, 'UniformOutput', false));
else
    AN = zeros(sn.nnodes, sn.nchains);
end
end
