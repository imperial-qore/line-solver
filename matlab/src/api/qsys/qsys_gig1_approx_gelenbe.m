function W=qsys_gig1_approx_gelenbe(lambda,mu,ca,cs)
% W=QSYS_GIG1_APPROX_GELENBE(LAMBDA,MU,CA,CS)

rho=lambda/mu;
W=(rho*ca^2+cs^2)/2/(1-rho)/lambda;
end
