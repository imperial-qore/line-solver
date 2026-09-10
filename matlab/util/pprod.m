function [n]=pprod(n,N)
% n=PPROD(N) % init
% n=PPROD(n,N) % next state
% Return a sequence of non-negative vectors less than a given vector
%
% n=pprod(N);
% while n>=0
%   ...
%   n=pprod(n,N);
% end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
    if nargin==1
        N=n;
        n=zeros(size(N));        
        return;
    end
    
    R=length(N);
    if sum(n==N)==R
        n=-1;
        return
    end

    s=R;
    while s>0 && n(s)==N(s)
        n(s)=0;
        s=s-1;
    end
    if s==0 
        %n=-1*ones(1,R);        
        return 
    end
    n(s)=n(s)+1;    
    return;
end