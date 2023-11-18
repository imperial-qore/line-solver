function [M,MU,MI,pi0]=cache_miss_asy(gamma,m,lambda)
% FPI method
[n,h]=size(gamma);
xi = cache_xi_fp(gamma,m);
MI=zeros(n,1);
for i=1:n
    MI(i) = sum(lambda(:,i,1))/(1+gamma(i,:)*xi(:));
end
M=sum(MI);
if nargin>2
    u=size(lambda,1);
    n=size(lambda,2);
    MU=zeros(u,1);
    for i=1:n
        pi0(i) = 1/(1+gamma(i,:)*xi(:));
        for v=1:u
            MU(v) = MU(v) + lambda(v,i,1)*pi0(i);
        end
    end
end
end