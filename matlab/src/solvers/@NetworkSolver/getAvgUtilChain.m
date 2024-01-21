function [UN] = getAvgUtilChain(self,U)
% [UN] = GETAVGUTILCHAIN(SELF,U)
% Return average utilization aggregated by chain
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
%if nargin == 1
%    [Q] = getAvgHandles(self);
%end
[UNclass] = getAvgUtil(self);

% compute average chain metrics
UN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        if ~isempty(UNclass)
            UN(i,c) = sum(UNclass(i,inchain)); %#ok<FNDSB>
        end
    end
end
end
