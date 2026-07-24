function [N]=oner(N,r)
% N=ONER(N,r)
% Decrement element in position of r of input vector
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
for s=r(:)'
    if s~=0
        N(s)=N(s)-1;
    end
end
end

