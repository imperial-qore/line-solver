function [W,rhohat]=qsys_gig1_approx_gelenbe(lambda,mu,ca,cs)
% [W,RHOHAT]=QSYS_GIG1_APPROX_GELENBE(LAMBDA,MU,CA,CS)
%
% Gelenbe's diffusion approximation with instantaneous-return boundary:
%   p(0) = 1-rho,  p(n) = rho*(1-rhat)*rhat^(n-1), n>=1
%   rhat = exp(-2*(1-rho)/(rho*ca^2+cs^2))
% hence E[N] = rho/(1-rhat) and the mean response time (time in system)
% W = E[N]/lambda = 1/(mu*(1-rhat)).
%
% Reference: Gelenbe, E. (1975). On approximate computer system models.
% Journal of the ACM 22(2), 261-269.

rho=lambda/mu;
rhat=exp(-2*(1-rho)/(rho*ca^2+cs^2));
W=1/(mu*(1-rhat));
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
