function [AMAP,AMAPS] = amap2_fit_gamma(M1, M2, M3, GAMMA) 
% Finds an AMAP(2) fitting the given characteristics.
% Input:
% - M1, M2, M3: moments of the inter-arrival times
% - GAMMA: auto-correlation decay rate of the inter-arrival times
% Output:
% - AMAP: the fitted AMAP(2)
% - AMAPS: all the fitted AMAP(2)

% if coefficient of variation is equal to 1, fit marked poisson
if abs(M2 - 2 * M1^2) < 1e-6
%   m3a_printf('Fitting AMAP(2): CV = 1, fitting a Poisson process\n');
   AMAP = {-1/M1, 1/M1};
   AMAPS = {AMAP};
   return;
end

% find all solutions for the parameters
AMAPS = amap2_fitall_gamma(M1, M2, M3, GAMMA);
for j = 1:length(AMAPS)
    AMAPS{j} = map_normalize(AMAPS{j});
end

% if no solutions is found, perform approximate fitting
if isempty(AMAPS)
    % find feasible characteristics
    [M2a, M3a, GAMMAa] = amap2_adjust_gamma(M1, M2, M3, GAMMA);
    % fit (should found at least one solution)
    AMAPS = amap2_fitall_gamma(M1, M2a, M3a, GAMMAa);
%     if isempty(AMAPS)
%         error('Fitting AMAP(2): feasibility could not be restored');
%     else
%         m3a_printf('Fitting AMAP(2): %d approximate solutions\n', length(AMAPS));
%         m3a_printf('Fitting AMAP(2): M2 = %f -> %f\n', M2, M2a);
%         m3a_printf('Fitting AMAP(2): M3 = %f -> %f\n', M3, M3a);
%         m3a_printf('Fitting AMAP(2): GAMMA = %f -> %f\n', GAMMA, GAMMAa);
%     end
else
%    m3a_printf('Fitting AMAP(2): %d exact solutions\n', length(AMAPS));
end

% if feasibility could not be restored (e.g. extreme moments from a
% near-null flow), fall back to a Poisson process matching the mean, as in
% the CV=1 branch above. This keeps the mean exact and yields a valid
% order-1 representation that downstream marked-fitting consumes.
if isempty(AMAPS)
    AMAP = {-1/M1, 1/M1};
    AMAPS = {AMAP};
    return;
end

AMAP = AMAPS{1};

end
