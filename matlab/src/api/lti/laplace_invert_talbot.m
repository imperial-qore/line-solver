function res = laplace_invert_talbot(F, t, n)
% RES = LAPLACE_INVERT_TALBOT(F, T, N)
%
% Invert a Laplace transform at T by Talbot's deformed contour. F is a handle
% called with a COMPLEX argument. N defaults to 32. Twin of the native Python
% api.lti.laplace_invert_talbot.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(n)
    n = 32;
end
if t <= 0
    line_error(mfilename, 'The Laplace inversion time point must be positive.');
end

alpha = talbot_get_alpha(n);
omega = talbot_get_omega(n, alpha);
res = 0.0;
for i = 1:n
    res = res + real(omega(i) * F(alpha(i)/t));
end
res = res / t;
end
