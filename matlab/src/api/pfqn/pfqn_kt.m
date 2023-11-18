function [G,lG,X,Q]=pfqn_kt(L,N,Z)
% Knessl-Tier asymptotic expansion
if isempty(L) || isempty(N) || sum(N)==0
    G=1;
    return;
end
if nargin<3
    Z=N*0;
end
[M,R]=size(L);
Ntot=sum(N);
for r=1:R
    beta(r)=N(r)/Ntot;
end
[X,Q]=pfqn_aql(L,N,Z);
delta=eye(R,R);
for i=1:R
    for j=1:R
        SK = 0;
        for k=1:M
            SK = SK + X(i)*X(j)*L(k,i)*L(k,j)/(1-sum(X(:)'*L(k,:)'))^2;
        end
        C(i,j) = delta(i,j)*beta(i) + (1/Ntot) * SK;
    end
end
Den=1;
for k=1:M
    Den=Den*max(1e-6,(1-sum(X(:)'*L(k,:)')));
end
%G=(2*pi).^(-R/2)/sqrt(Ntot^R*det(C))*exp(-Ntot*beta*log(X)')/Den;
lG = log((2*pi).^(-R/2)/sqrt(Ntot^R*det(C))) + (-Ntot*beta*log(X)') - log(Den);
G=exp(lG);
end