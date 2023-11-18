function [Gn,lGn]=pfqn_panacea(L,K,Z)
% [GN,LGN]=PFQN_PANACEA(L,K,Z)

% K = population vector
[q,p]=size(L);
if nargin==2 || isempty(Z)
    Z=K*0+1e-8;
end
if isempty(L) | sum(L,1)==zeros(1,p)
    lGn = - sum(factln(K)) + sum(K.*log(sum(Z,1)));
    Gn=exp(lGn);
    return
end
r = L./repmat(Z,q,1);
N = max(1./r(:));
beta = K/N;
gamma = r * N;
alpha = 1-K*r';
gammatilde = gamma ./ repmat(alpha',1,p);
if min(alpha)<0
    %    line_warning(mfilename,'Model is not in normal usage');
    Gn=NaN;
    lGn=NaN;
    return
end

A0 = 1;
A1 = 0;
for j=1:p
    m = zeros(1,p); m(j)=2;
    A1 = A1 -beta(j) * pfqn_ca(gammatilde,m);
end

A2 = 0;
for j=1:p
    m = zeros(1,p); m(j)=3;
    A2 = A2 + 2 * beta(j) * pfqn_ca(gammatilde,m);
    m = zeros(1,p); m(j)=4;
    A2 = A2 + 3 * beta(j)^2 * pfqn_ca(gammatilde,m);
    for k=setdiff(1:p,j)
        m = zeros(1,p); m(j)=2; m(k)=2;
        A2 = A2 + 0.5 * beta(j) * beta(k) * pfqn_ca(gammatilde,m);
    end
end

if 0
    A3 = 0;
    for j=1:p
        m = zeros(1,p); m(j)=4;
        A3 = A3 - 6 * beta(j) * pfqn_ca(gammatilde,m);
        m = zeros(1,p); m(j)=5;
        A3 = A3 - 20 * beta(j)^2 * pfqn_ca(gammatilde,m);
        m = zeros(1,p); m(j)=6;
        A3 = A3 - 15 * beta(j)^3 * pfqn_ca(gammatilde,m);
        for k=setdiff(1:p,j)
            m = zeros(1,p); m(j)=4; m(k)=2;
            A3 = A3 - 2 * beta(j) * beta(k) * pfqn_ca(gammatilde,m);
            m = zeros(1,p); m(j)=2; m(k)=3;
            A3 = A3 - 3 * beta(j)^2 * beta(k) * pfqn_ca(gammatilde,m);
            for l=setdiff(1:p,[j,k])
                m = zeros(1,p); m(j)=2; m(k)=2; m(l)=2;
                A3 = A3 - (1/6) * beta(j) * beta(k) * beta(l) * pfqn_ca(gammatilde,m);
            end
        end
    end
end
I = [A0, A1/N, A2/N^2];%, A3/N^3*0];

lGn = - sum(factln(K)) + sum(K.*log(sum(Z,1))) + log(sum(I)) - sum(log(alpha));
Gn = exp(lGn);
if ~isfinite(lGn)
    Gn=NaN;
    lGn=NaN;
end
end
