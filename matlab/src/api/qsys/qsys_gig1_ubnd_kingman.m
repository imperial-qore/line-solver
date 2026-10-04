function [W,rhohat]=qsys_gig1_ubnd_kingman(lambda,mu,ca,cs)
% [W,RHOHAT]=QSYS_GIG1_UBND_KINGMAN(LAMBDA,MU,CA,CS)
%
% Kingman's upper bound on the mean waiting time of a G/G/1 queue:
%   Wq <= lambda*(sa^2+ss^2)/(2*(1-rho)),  sa^2=ca^2/lambda^2, ss^2=cs^2/mu^2
% The returned W adds the mean service time, so it upper-bounds the mean
% response time (time in system).
%
% Reference: Kingman, J.F.C. (1962). Some inequalities for the queue GI/G/1.
% Biometrika 49(3/4), 315-324.

rho=lambda/mu;
Wq = lambda*(ca^2/lambda^2 + cs^2/mu^2)/(2*(1-rho));
W = Wq + 1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
