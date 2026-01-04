function MAP = amap2_fit_gamma_map(map) 
% Performs approximate fitting of a given MAP, yielding a second-order
% MAP in canonical form.
%
% Input
% - map: the MAP (of arbitrary order) to fit
% Output
% - MAP: fitted second-order MAP

M1 = map_mean(map);
M2 = map_moment(map,2);
M3 = map_moment(map,3);
GAMMA = map_gamma(map);

MAP = amap2_fit_gamma(M1, M2, M3, GAMMA);

end