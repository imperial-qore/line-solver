function APH = aph2_fit_map(map) 
% Performs approximate fitting of a MAP, yielding a second-order
% APH in canonical form.
% Input
% - map: the MAP of arbitrary order to fit
% Output
% - APH: fitted second-order phase-type

M1 = map_mean(map);
M2 = map_moment(map,2);
M3 = map_moment(map,3);

APH = aph2_fit(M1, M2, M3);

end