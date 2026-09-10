function res = laplace_invert_gaver_stehfest(F, t, n)
% RES = LAPLACE_INVERT_GAVER_STEHFEST(F, T, N)
%
% Invert a Laplace transform at T by Gaver-Stehfest. This method samples the
% REAL axis only, so F may be a real-argument handle. N defaults to 12 and is
% rounded DOWN to even. Twin of the native Python
% api.lti.laplace_invert_gaver_stehfest.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(n)
    n = 12;
end
if mod(n,2) == 1
    n = n - 1;
end
if t <= 0
    line_error(mfilename, 'The Laplace inversion time point must be positive.');
end

alpha = gaver_stehfest_get_alpha(n);
omega = gaver_stehfest_get_omega(n);
res = 0.0;
for i = 1:n
    res = res + omega(i) * F(alpha(i)/t);
end
res = res / t;
end
