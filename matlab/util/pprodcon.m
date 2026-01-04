function [n]=pprodcon(n,lb,ub)
% n=PPRODCON(N,lb,ub) - init
% n=PPRODCON(n,N)
% Runs pprod with constraints on vector range 
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% initialization
if nargin==2
    lb=n;
    n=lb;
    return;
end

% stop condition
R=length(ub);
if sum(n==ub)==R
    n=-1;
    return
end

s=R;
while s>0 && n(s)==ub(s)
    n(s)=lb(s);
    s=s-1;
end
if s==0
    %n=-1*ones(1,R);
    return
end
n(s)=n(s)+1;
return;
end