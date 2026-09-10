function res = euler_get_omega(n)
% RES = EULER_GET_OMEGA(N)
%
% Euler (Abate-Whitt) weights omega_i = 10^((n-1)/6) (-1)^(i-1) eta_i. Twin of
% the native Python api.lti.euler_get_omega.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

eta = euler_get_eta(n);
res = (10.0^((n-1)/6.0)) * ((-1.0).^(0:n-1)) .* eta;
end
