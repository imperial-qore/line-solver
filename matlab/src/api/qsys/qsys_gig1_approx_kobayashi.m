function [W,rhohat]=qsys_gig1_approx_kobayashi(lambda,mu,ca,cs)
% W=QSYS_GIG1_APPROX_KOBAYASHI(LAMBDA,MU,CA,CS)

rho=lambda/mu;
rhohat=exp(-2*(1-rho)/(rho*(ca^2+cs^2/rho)));
W=rhohat/(1-rhohat)/lambda;
end
