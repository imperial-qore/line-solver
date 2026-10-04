function [W,rhohat]=qsys_gig1_approx_myskja2(lambda,mu,ca,cs,q0,qa)
% [W,RHOHAT]=QSYS_GIG1_APPROX_MYSKJA2(LAMBDA,MU,CA,CS,Q0,QA)
%
% Myskja's enhanced third-moment approximation of the mean response time
% (time in system). For ca=1 the interpolation parameter theta is a 0/0
% form, so the exact M/G/1 result is returned instead.
%
% qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
% q0 = lowest value of the relative third moment for a given mean and SCV

if abs(ca^2-1) < 1e-8
    % M/G/1 case: exact (also the interpolation anchor of the method)
    [W,rhohat] = qsys_mg1(lambda,mu,cs);
    return
end
ra = (1+ca^2)/2;
rs = (1+cs^2)/2;
rho=lambda/mu;
theta=(rho*(qa-ra)-(qa-ra^2))/(2*rho*(ra-1));
d=(1+1/ra)*(1-rs)*(1-(q0/qa)^3)*(1-rho^3);
D = (rs-theta)^2+(2*rs-1+d)*(ra-1);
D = max(D,0); % guard small negative values due to round-off
W=(rho/(1-rho))/lambda*(rs+(1/rho)*(sqrt(D)-(rs-theta)));
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
