function [out1,out2,out3,out4]=sprod(par1,par2,par3)
% [s,n,S,D]=SPROD(M,N) % init
% [s,n]=SPROD(s,S,D) % next state
% Sequence of non-negative matrices with constant row sums
% 
% Example: 
% [s,n,S,D]=sprod(2,5); % m=2 elements summing to N=5
%
% while s>=0
% ... 
% [s,n]=sprod(s,S,D) % generate next state
% end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin==2 % this is init
    M=par1; N=par2;
    R=length(N);
    S=0;
    for r=1:R
        D{r}=multichoose(M,N(r));
        S(r)=size(D{r},1);
    end
    S=S-1;
    s=pprod(S);
    for r=1:R
        n(:,r) = D{r}(1+s(r),:)';
    end
else
    s=par1; S=par2; D=par3;
    s=pprod(s,S);
    M=size(D{1},2);
    R=size(D,2);
    n=zeros(M,R);
    if s==-1
        n=n-1;
        out1=s;
        out2=n;
        out3=S;
        out4=D;
        return
    end
    for r=1:R
        n(:,r) = D{r}(1+s(r),:)';
    end
end
out1=s;
out2=n;
out3=S;
out4=D;
end