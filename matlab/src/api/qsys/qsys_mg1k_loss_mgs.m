function [lossprob_mg1k,rho]=qsys_mg1k_loss_mgs(lambda,mu,mu_scv,K)
% J. MacGregor Smith - Optimal Design and Performance Modelling of M/G/1/K Queueing Systems
rho = lambda/mu;
s = sqrt(mu_scv);
sqrt_rho = sqrt(rho);
lossprob_mg1k_num = rho^((sqrt_rho*s^2-sqrt_rho+2*K)/(2+sqrt_rho*s^2-sqrt_rho))*(rho-1);
lossprob_mg1k_den = (rho^(2*(1+sqrt_rho*s^2-sqrt_rho+K)/(2+sqrt_rho*s^2-sqrt_rho))-1);
lossprob_mg1k = lossprob_mg1k_num  / lossprob_mg1k_den;
end
