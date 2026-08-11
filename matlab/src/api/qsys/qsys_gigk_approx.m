function [W,rhohat]=qsys_gigk_approx(lambda,mu,ca,cs,k)
% W=QSYS_GIGK_APPROX(LAMBDA,MU,CA,CS,K)
rho = lambda / (mu*k);
if rho>0.7
    alpha = (rho^k+rho)/2;
else
    alpha = rho^((k+1)/2);
end
W = (alpha/mu)*(1/(1-rho))*(ca^2+cs^2)/(2*k) + 1/mu; % includes service time
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
