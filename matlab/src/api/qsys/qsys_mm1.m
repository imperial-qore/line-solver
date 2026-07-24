function [W,rho]=qsys_mm1(lambda,mu)
% W=QSYS_MM1(LAMBDA,MU)
rho=lambda/mu;
W=rho/(1-rho) / lambda;
end
