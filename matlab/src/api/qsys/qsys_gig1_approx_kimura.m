function W=qsys_gig1_approx_kimura(sigma,mu,ca,cs)
% W=QSYS_GIG1_APPROX_KIMURA(SIGMA,MU,CA,CS)

W=sigma*(ca^2+cs^2)/mu/(1-sigma)/(1+ca^2);
end
