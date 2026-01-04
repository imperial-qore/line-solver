function W=qsys_gig1_approx_myskja2(lambda,mu,ca,cs,q0,qa)
% W=QSYS_GIG1_APPROX_MYSKJA2(LAMBDA,MU,CA,CS,Q0,QA)

% qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
% q0 = lowest value of the relative third moment for a given mean and SCV
ra = (1+ca^2)/2;
rs = (1+cs^2)/2;
rho=lambda/mu;
theta=(rho*(qa-ra)-(qa-ra^2))/(2*rho*(ra-1));
d=(1+1/ra)*(1-rs)*(1-(q0/qa)^3)*(1-rho^3);
D = (rs-theta)^2+(2*rs-1+d)*(ra-1);
W=(rho/(1-rho))/lambda*(rs+(1/rho)*(sqrt(D)-(rs-theta)));
end
