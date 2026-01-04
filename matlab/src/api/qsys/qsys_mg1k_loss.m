function [lossprob_mg1k,rho]=qsys_mg1k_loss(lambda,svc_density,K)
% Niu-Cooper, Transform-Free Analysis of M/G/1/K and Related Queues, Mathematics of Operations Research
% Vol. 18, No. 2 (May, 1993), pp. 486-510 (25 pages)
% TEST:
% mu=3; lambda=2; K=5; [lossprob_mg1k,rho]=qsys_mg1k_loss(lambda,@(t)mu.*exp(-mu.*t),K)
%
%sigma is the embedded probability seen immediately after the kth *service-start* epoch
tmax = 1e4/lambda;
mu = 1/integral(@(t) t.*svc_density(t),0,tmax);
a = zeros(1,K-1);
for j=0:1e3
    factj = factorial(j);
    a(1+j)=integral(@(t)exp(-lambda*t).*((lambda*t).^j).*svc_density(t),0,tmax)/factj;
    if a(1+j)<1e-12
        break
    end
end
% We model the number of customers waiting in the queue immediately
% after the kth *service-start* epoch
P = zeros(K-1);
if false
    P(1,1) = a(1+0)+a(1+1);
    P(2,1) = a(1+0);
    for j=1:(K-3)
        for i=0:(j+1)
            P(1+i,1+j) = a(1+(j+1-i));
        end
    end
    for i=0:(K-2)
        P(1+i,1+(K-2)) = sum(a(1+((K-2+1-i)):end));
    end
else
    % traditional embedded process at departure
    for i=0:(K-2)
        P(1,1+i)=a(1+i);
        P(2,1+i)=a(1+i);
    end
    P(1,K-1)=1-sum(P(1,:));
    P(2,K-1)=1-sum(P(2,:));
    for j=3:K-1
        for i=j-1:(K-2)
            P(j,i)=a(1+i-j+1);
        end
        P(j,K-1)=1-sum(P(j,1:K-2));
    end
end
P=dtmc_makestochastic(P);
sigma = dtmc_solve(P);
rho = lambda / mu;
lossprob_mg1k = 1-1/(sigma(1+0)*a(1+0)+rho);
end
