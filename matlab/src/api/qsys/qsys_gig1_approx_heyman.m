function [W,rhohat]=qsys_gig1_approx_heyman(lambda,mu,ca,cs)
% W=QSYS_GIG1_APPROX_HEYMAN(LAMBDA,MU,CA,CS)

rho=lambda/mu;
W=rho/(1-rho)/mu*(ca^2+cs^2)/2+1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
