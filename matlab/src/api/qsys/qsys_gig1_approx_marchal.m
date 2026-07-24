function [W,rhohat]=qsys_gig1_approx_marchal(lambda,mu,ca,cs)
% W=QSYS_GIG1_APPROX_MARCHAL(LAMBDA,MU,CA,CS)

rho=lambda/mu;
Wmm1 = rho/(1-rho);
W = Wmm1*(1+cs^2)/2/mu*(ca+rho^2*cs^2)/(1+rho^2*cs^2)+1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
