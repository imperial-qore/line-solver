function [AN] = getAvgArvRChain(self,A)
% [AN] = GETAVGARVRCHAIN(SELF,A)
% Return average arrival rates aggregated by chain
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
if nargin == 1
    [~,~,~,~,A] = getAvgHandles(self);
end
[ANclass] = getAvgArvR(self);

% compute average chain metrics
AN = zeros(sn.nstations, sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for i=1:sn.nstations
        if ~isempty(ANclass)
            AN(i,c) = sum(ANclass(i,inchain));
        end
    end
end
end
