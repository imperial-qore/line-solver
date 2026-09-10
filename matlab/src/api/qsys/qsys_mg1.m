function [W,rhohat]=qsys_mg1(lambda,mu,cs)
% W=QSYS_MG1(LAMBDA,MU,CS)

rho=lambda/mu;
Q = rho + rho^2/(2*(1-rho)) + lambda^2*cs^2/mu^2/(2*(1-rho));
W = Q/lambda;
rhohat = Q/(1+Q); % so that M/M/1 formulas still hold
end
