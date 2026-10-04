function [TN] = getAvgTputChain(self,T)
% [TN] = GETAVGTPUTCHAIN(SELF,T)

% Return average throughputs aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
TNclass = getAvgTput(self);

% compute average chain metrics
if ~isempty(TNclass)
    TN = cell2mat(cellfun(@(chain) sum(TNclass(:, chain), 2), sn.inchain, 'UniformOutput', false));
else
    TN = zeros(sn.nstations, sn.nchains);
end

end
