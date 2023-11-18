function [lossprob_mm1k,rho]=qsys_mm1k_loss(lambda,mu,K)
% Niu-Cooper, Transform-Free Analysis of M/G/1/K and Related Queues, Mathematics of Operations Research
% Vol. 18, No. 2 (May, 1993), pp. 486-510 (25 pages)
% TEST:
% mu=3; lambda=2; K=5; [sigma,rho,lossprob]=qsys_mg1k_loss(lambda,@(t)mu.*exp(-mu.*t),K)
rho = lambda/mu;
lossprob_mm1k = (1-rho)/(1-rho^(K+1))*rho^K;
end
