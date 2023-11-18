function [QN] = getAvgQLenChain(self,Q)
% [QN] = GETAVGQLENCHAIN(SELF,Q)

% Return average queue-lengths aggregated by chain
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
%if nargin == 1
%    [Q] = getAvgHandles(self);
%end
[QNclass] = getAvgQLen(self);

% compute average chain metrics
QN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        if ~isempty(QNclass)
            QN(i,c) = sum(QNclass(i,inchain));
        end
    end
end
end
