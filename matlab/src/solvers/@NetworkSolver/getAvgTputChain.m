function [TN] = getAvgTputChain(self,T)
% [TN] = GETAVGTPUTCHAIN(SELF,T)

% Return average throughputs aggregated by chain
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
%if nargin == 1
%    [Q] = getAvgHandles(self);
%end
[TNclass] = getAvgTput(self);

% compute average chain metrics
TN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        if ~isempty(TNclass)
            TN(i,c) = sum(TNclass(i,inchain)); %#ok<FNDSB>
        end
    end
end
end
