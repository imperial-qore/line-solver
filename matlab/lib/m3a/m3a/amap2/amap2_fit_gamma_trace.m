function mmap = amap2_fit_gamma_trace(T) 
% Performs approximate fitting of a given trace, yielding a second-order
% MAP in canonical form.
% Input
% - T: the inter-arrival times
% Output
% - mmap: fitted MAMAP[m]

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);
GAMMA = trace_gamma(T);

mmap = amap2_fit_gamma(M1, M2, M3, GAMMA);

end