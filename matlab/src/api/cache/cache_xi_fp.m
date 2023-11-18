function [xi,pi0,pij,it] = cache_xi_fp(gamma,m,xi)
[n,h]=size(gamma);
tol=1e-14;
pi0=ones(1,n)/(h+1);
pij=zeros(n,h);
xi=zeros(1,h);
if nargin<3
    for l=1:h
        xi(l) = m(l)/mean(gamma(:,l))/(n+sum(m)-1);
    end
end
for it=1:1e4
    pi0_1=pi0;
    xi = m ./ (pi0_1*gamma);
    pij = abs(gamma .* repmat(xi,n,1)) ./ abs(1+gamma*xi');
    pi0=max(tol,1-sum(pij,2)');
    DELTA=norm(abs(1-pi0(:)./pi0_1(:)),1);
    if DELTA<tol
        xi(xi<0)=tol;
        return
    end
end
xi(xi<0)=tol;
end
