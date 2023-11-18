function APH = aph2_fit_trace(T) 
% Performs approximate fitting of a given trace, yielding a second-order
% APH in canonical form.
% Input
% - T: the inter-arrival times
% Output
% - APH: fitted second-order phase-type

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);

APH = aph2_fit(M1, M2, M3);

end