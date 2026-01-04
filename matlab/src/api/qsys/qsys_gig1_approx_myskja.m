function W=qsys_gig1_approx_myskja(lambda,mu,ca,cs,q0,qa)
% W=QSYS_GIG1_APPROX_MYSKJA(LAMBDA,MU,CA,CS,Q0,QA)

% qa = third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
% q0 = lowest value of the relative third moment for a given mean and SCV
rho=lambda/mu;
W=rho/(2*mu*(1-rho))*((1+cs^2)+(q0/q)^(1/rho-rho)*(1/rho)*(ca^2-1));
end
