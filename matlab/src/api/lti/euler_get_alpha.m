function res = euler_get_alpha(n)
% RES = EULER_GET_ALPHA(N)
%
% Euler (Abate-Whitt) nodes alpha_i = (n-1)log(10)/6 + i*pi*1i, i = 0..N-1.
% Twin of the native Python api.lti.euler_get_alpha.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

res = complex((n-1) * log(10.0) / 6, pi * (0:n-1));
end
