function [W,rhohat]=qsys_gig1_approx_myskja(lambda,mu,ca,cs,q0,qa)
% [W,RHOHAT]=QSYS_GIG1_APPROX_MYSKJA(LAMBDA,MU,CA,CS,Q0,QA)
%
% Myskja's third-moment approximation of the mean waiting time
%   Wq = rho/(2*mu*(1-rho))*((1+cs^2)+(q0/qa)^(1/rho-rho)*(1/rho)*(ca^2-1))
% exact for M/G/1 (ca=1). The returned W adds the mean service time
% (response time, time in system).
%
% qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
% q0 = lowest value of the relative third moment for a given mean and SCV

rho=lambda/mu;
Wq=rho/(2*mu*(1-rho))*((1+cs^2)+(q0/qa)^(1/rho-rho)*(1/rho)*(ca^2-1));
W=Wq+1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
