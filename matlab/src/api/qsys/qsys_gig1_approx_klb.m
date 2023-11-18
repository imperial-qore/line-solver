function [W,rhohat]=qsys_gig1_approx_klb(lambda,mu,ca,cs)
% [W,rhohat]=QSYS_GIG1_APPROX_KLB(LAMBDA,MU,CA,CS)

% kramer-langenbach-belz formula
rho=lambda/mu;
if ca<=1
    g=exp(-2*(1-rho)*(1-ca^2)^2/(3*rho*(ca^2+cs^2)));
else
    g=exp(-(1-rho)*(ca^2-1)/(ca^2+4*cs^2));
end
W=1/mu*((rho/(1-rho))*((cs^2+ca^2)/2)*g+1);
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
