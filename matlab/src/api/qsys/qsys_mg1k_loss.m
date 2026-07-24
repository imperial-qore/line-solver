function [lossprob_mg1k,rho]=qsys_mg1k_loss(lambda,svc_density,K)
% [LOSSPROB,RHO]=QSYS_MG1K_LOSS(LAMBDA,SVC_DENSITY,K)
%
% Exact M/G/1/K loss probability via the Markov chain embedded at
% service-start epochs (transform-free analysis in the spirit of
% Niu-Cooper).
%
% State: number of customers waiting in the queue immediately after a
% service start, q in {0,...,K-2} (capacity K includes the job in service;
% just after a departure at most K-1 jobs remain, one of which enters
% service). With a_j = P(j Poisson arrivals during a service time):
%   q=0 : if no arrival occurs during the service the system empties and
%         the next service starts with the next arrival (q'=0), so both
%         a_0 and a_1 lead to q'=0 and j>=2 arrivals lead to q'=j-1;
%   q>=1: q' = q-1+j, with arrivals beyond the free capacity lost
%         (aggregated in the last column).
% The loss probability follows from the renewal-reward argument
%   E[cycle] = E[S] + sigma_0*a_0/lambda,  lambda_eff = 1/E[cycle],
%   P_loss = 1 - lambda_eff/lambda = 1 - 1/(rho + sigma_0*a_0)
% where sigma is the stationary distribution at service-start epochs.
%
% TEST:
% mu=3; lambda=2; K=5; [lossprob_mg1k,rho]=qsys_mg1k_loss(lambda,@(t)mu.*exp(-mu.*t),K)
% (exact M/M/1/5 value: 0.04812030)
%
% Reference: Niu, Cooper. Transform-Free Analysis of M/G/1/K and Related
% Queues. Mathematics of Operations Research 18(2), 1993, 486-510.

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
% Embedded chain at service-start epochs, states q=0..K-2 (row/col q+1)
P = zeros(K-1);
% row 1 (q=0): idle period after an empty departure epoch
P(1,1) = a(1+0)+a(1+1);
for i=1:(K-3)
    P(1,1+i) = a(1+i+1);
end
P(1,K-1) = 1-sum(P(1,1:K-2));
% row 2 (q=1): q' = number of arrivals during the service (capped)
if K >= 3
    for i=0:(K-3)
        P(2,1+i)=a(1+i);
    end
    P(2,K-1)=1-sum(P(2,1:K-2));
end
% rows j>=3 (q=j-1): q' = q-1+arrivals (capped)
for j=3:K-1
    for i=j-1:(K-2)
        P(j,i)=a(1+i-j+1);
    end
    P(j,K-1)=1-sum(P(j,1:K-2));
end
P=dtmc_makestochastic(P);
sigma = dtmc_solve(P);
rho = lambda / mu;
lossprob_mg1k = 1-1/(sigma(1+0)*a(1+0)+rho);
end
