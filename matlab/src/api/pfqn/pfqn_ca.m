%{
%{
 % @file pfqn_ca.m
 % @brief Convolution Algorithm for exact normalizing constant computation.
%}
%}

function [Gn,lGn]=pfqn_ca(L,N,Z)
%{
%{
 % @brief Convolution Algorithm for exact normalizing constant computation.
 % @fn pfqn_ca(L, N, Z)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @return Gn Normalizing constant.
 % @return lGn Logarithm of the normalizing constant.
%}
%}
[M,R]=size(L);
if M==0
    lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1)));
    Gn = exp(lGn);
    return
end

if min(N)<0
    Gn=0;
    lGn=-Inf;
    return;
end

if sum(N)==0
    Gn=1;
    lGn=0;
    return;
end


if nargin<3 || isempty(Z)
    Z=zeros(1,R);
end

G = ones(M+1,prod(N+1)); % stores G across recursion
n = pprod(N);
while sum(n)~=-1
    idxn = hashpop(n,N);
    G(1,idxn) = Fz(Z,n);
    for m=2:M+1
        G(m,idxn) = G(m-1,idxn); % norm constant with m-1 queues
        for r=1:R
            if n(r)>=1
                n(r) = n(r)-1;
                idxn_1r = hashpop(n,N);
                n(r) = n(r)+1;
                G(m,idxn) = G(m,idxn) + L(m-1,r)*G(m,idxn_1r);
            end
        end
    end
    n=pprod(n,N);
end
Gn=G(M+1,end);
lGn=log(Gn);
end

function idx=hashpop(n,N,R,prods)
% IDX=HASHPOP(N,N,R,PRODS)

% hash a population vector in n: 0<=n<=N
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

function [n]=pprod(n,N)
% [N]=PPROD(N,N)

% sequentially generate all vectors n: 0<=n<=N
% n=pprod(N) - init
% n=pprod(n,N) - next state
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

function f=Fz(Z,n)
% F=FZ(Z,N)

R=length(n);
if sum(n)==0
    f=1;
    return
end
f=0;
for r=1:R
    if Z(r)>0
        f=f+log(Z(r))*n(r);
        f=f-gammaln(1+n(r));
    elseif n(r)>0
        f = 0;
        return
    end
end
f=exp(f);
end
