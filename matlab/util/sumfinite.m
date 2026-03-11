function s = sumfinite(v, dim)
% S = SUMFINITE(V, DIM)
% Sum the finite values in vector V alonside dimensions DIM
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

v(~isfinite(v)) = 0;
if nargin>1
    if dim == 1
        s = sum(v, 1);
    elseif dim == 2
        s = sum(v, 2);
    else
        s = sum(v);
    end
else
    s = sum(v);
end