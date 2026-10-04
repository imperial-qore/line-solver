function idx=hashpop(n,N,R,prods)
% idx=HASHPOP(n,N)
% idx=HASHPOP(n,N,R,prods) where prods(r)=prod(N(1:r-1)+1) (faster)
% Hash a vector n on an integer lattice defined by vector N
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
idx=1;
if nargin==2
    R=length(N);
    for r=1:R
        idx= idx + prod(N(1:r-1)+1)*n(r);
    end
    return
else
    for r=1:R
        idx= idx + prods(r)*n(r);
    end
end
end