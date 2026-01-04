function W=qsys_gm1(sigma,mu)
% W=QSYS_GM1(SIGMA,MU)

% sigma = Load at arrival instants (Laplace transform of the inter-arrival times)
W=1/(1-sigma)/mu;
end
