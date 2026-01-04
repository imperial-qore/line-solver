%{
%{
 % @file pfqn_kt.m
 % @brief Knessl-Tier asymptotic expansion for normalizing constant.
%}
%}

%{
%{
 % @brief Knessl-Tier asymptotic expansion for normalizing constant.
 % @fn pfqn_kt(L, N, Z)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector (default: zeros).
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
 % @return X System throughput.
 % @return Q Mean queue lengths.
%}
%}
function [G,lG,X,Q]=pfqn_kt(L,N,Z)
% Knessl-Tier asymptotic expansion
if isempty(L) || isempty(N) || sum(N)==0
    G = 1;
    return;
end
if nargin<3
    Z = N*0;
end
[Morig,Rorig] = size(L);
% fix self-looping customers as they would yield Uk=1
%slcs = 0;
slcdemandfactor = 0;
if Rorig>1
    isslc = false(1,Rorig);
    for r=1:Rorig
        if nnz(L(:,r))==1
            if Z(r)==0
                ist = find(L(:,r)>0);
                L = [L; repmat(L(ist,:),N(r),1)];
                isslc(r) = true;
                slcdemandfactor = N(r)*log(L(ist,r));
            end
            %slcs = slcs + N(r);
        end
    end
    L(:,isslc)=[];
    Z(:,isslc)=[];
    N(:,isslc)=[];
end
[M,R] = size(L);
Ntot = sum(N);
beta = zeros(1,R);
for r=1:R
    beta(r) = N(r)/Ntot;
end
if Ntot <= 4
    [X,Q] = pfqn_bs(L,N,Z);
else
    [X,Q] = pfqn_aql(L,N,Z);
end
delta = eye(R,R);
C = zeros(R);
for s=1:R
    for r=1:R
        SK = 0;
        for k=1:M
            SK = SK + X(s)*X(r)*L(k,s)*L(k,r)/max(GlobalConstants.FineTol,1-sum(X(:)'*L(k,:)'))^2;
        end
        C(s,r) = delta(s,r)*beta(s) + (1/Ntot) * SK;
    end
end
Den = 1;
for k=1:M
    Den = Den * max(GlobalConstants.FineTol, 1-sum(X(:)'*L(k,:)'));
end
%G=(2*pi).^(-R/2)/sqrt(Ntot^R*det(C))*exp(-Ntot*beta*log(X)')/Den;
lG = log((2*pi).^(-R/2)/sqrt(Ntot^R*det(C))) + (-Ntot*beta*log(X)') - log(Den) + slcdemandfactor;
G = exp(lG);
end