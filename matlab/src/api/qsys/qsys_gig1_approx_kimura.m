function [W,rhohat]=qsys_gig1_approx_kimura(lambda,mu,ca,cs)
% [W,RHOHAT]=QSYS_GIG1_APPROX_KIMURA(LAMBDA,MU,CA,CS)
%
% Kimura's diffusion-interpolation approximation of the mean waiting time
%   Wq = rho*(ca^2+cs^2)/(mu*(1-rho)*(1+ca^2))
% exact for M/M/1 and M/G/1. The returned W adds the mean service time
% (response time, time in system).
%
% Reference: Kimura, T. (1986). A two-moment approximation for the mean
% waiting time in the GI/G/s queue. Management Science 32(6), 751-763.

rho=lambda/mu;
Wq=rho*(ca^2+cs^2)/mu/(1-rho)/(1+ca^2);
W=Wq+1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
