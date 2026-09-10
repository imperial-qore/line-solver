function [UN] = getAvgNodeUtilChain(self,U)
% [UN] = GETAVGNODEUTILCHAIN(SELF,U)
% Return average utilization aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
[~, UNclass] = getAvgNode(self);

% compute average chain metrics
if ~isempty(UNclass)
    UN = cell2mat(cellfun(@(chain) sum(UNclass(:, chain), 2), sn.inchain, 'UniformOutput', false));
else
    UN = zeros(sn.nnodes, sn.nchains);
end
end
