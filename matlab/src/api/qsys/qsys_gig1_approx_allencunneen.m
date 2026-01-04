function [W,rhohat]=qsys_gig1_approx_allencunneen(lambda,mu,ca,cs)
% W=QSYS_GIG1_APPROX_ALLENCUNNEEN(LAMBDA,MU,CA,CS)

rho=lambda/mu;
W=(rho/(1-rho))/mu*((cs^2+ca^2)/2) + 1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
