function [M,MU,MI,pi0,lE]=cache_miss_rayint(gamma,m,lambda)
% M: global miss rate
% MU: per-user miss rate
% MI: per-item miss rate
ma=m; ma(1)=ma(1)+1;
[~,lE,xi] = cache_rayint(gamma,m);
[~,lEa,~] = cache_rayint(gamma,ma);
M =  exp(lEa - lE);
if nargin>2
    u=size(lambda,1);
    n=size(lambda,2);
    pi0=zeros(1,n);
    % compute MU
    MU=zeros(u,1);
    if nargout>1
        for k=1:n
            if sum(gamma(k,:))>0
                [~,lE1(k)]=cache_rayint(gamma(setdiff(1:n,k),:),m,xi);                
                pi0(k) = exp(lE1(k) - lE);
                if pi0(k)>1 || pi0(k)<0 % try to recompute xi
                    [~,lE1(k)]=cache_rayint(gamma(setdiff(1:n,k),:),m);                
                    pi0(k) = exp(lE1(k) - lE);
                end
                for v=1:u
                    MU(v) = MU(v) + (lambda(v,k,1))*pi0(k);
                end
            end
        end
    end
    % compute MI
    if nargout>2
        MI=zeros(n,1);
        for k=1:n
            if sum(gamma(k,:))>0
                MI(k) = MI(k) + sum(lambda(:,k,1))*exp(lE1(k)-lE);
            else
                MI(k)=0;
            end
        end
    end
end
end