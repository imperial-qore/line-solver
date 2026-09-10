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
if nargin<3 || isempty(Z)
    Z=zeros(1,R);
end
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

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
Nt = sum(N);
% Each class independently takes whichever station -- or the delay -- gives it
% its largest factor. The mixed state so named has term at least the product of
% those factors, because a station holding several classes carries a multinomial
% coefficient of at least one, so this is still a LOWER bound on log G. It
% dominates the per-configuration maximum it replaces, which asked ONE station
% (or the delay) to hold every class at once and so dropped the delay entirely
% as soon as a single class had no think time. That collapse is what made the
% scaling scale UP: on L=[1e-9,1], N=[99,1], Z=[1,0] the old estimate was the
% all-at-the-queue -2051.6 against a true log G of -359.1, giving kscale=-30,
% and Z/2^-30 = 1.07e9 overflowed the delay column Z^n/n! at n=[40,0].
lGest = 0;
for r = 1:R
    if N(r) <= 0, continue; end
    best = -Inf;
    for i = 1:M
        if L(i,r) > 0
            best = max(best, N(r)*log(L(i,r)));
        end
    end
    if sum(Z(:,r)) > 0
        best = max(best, N(r)*log(sum(Z(:,r))) - factln(N(r)));
    end
    if ~isfinite(best)
        % no station and no delay can hold class r, so G(N) is exactly zero
        lGest = -Inf; break
    end
    lGest = lGest + best;
end
if ~isfinite(lGest)
    kscale = 0;
else
    kscale = round(lGest/(Nt*log(2)));
end
cscale = pow2(kscale);
L = L / cscale;
Z = Z / cscale;

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
% Undo the scaling in log space: log G = log G_scaled + sum(N) log c. lGn is
% therefore finite whenever log G itself is, even though Gn may legitimately
% overflow to Inf (the true constant really is outside double range).
lGn = log(G(M+1,end)) + sum(N)*kscale*log(2);
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
Gn = pow2(G(M+1,end), sum(N)*kscale);
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
