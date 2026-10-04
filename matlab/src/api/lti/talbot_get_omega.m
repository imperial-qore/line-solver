function arr = talbot_get_omega(n, alpha)
% ARR = TALBOT_GET_OMEGA(N, ALPHA)
%
% Talbot contour weights. Twin of the native Python api.lti.talbot_get_omega.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(alpha)
    alpha = talbot_get_alpha(n);
end
arr = zeros(1, n);
arr(1) = exp(alpha(1)) / 5.0;
for i = 2:n
    theta = (i-1) * pi / n;
    cotTheta = 1.0 / tan(theta);
    multiplier = complex(1.0, theta * (1 + cotTheta^2) - cotTheta);
    arr(i) = 2 * exp(alpha(i)) / 5.0 * multiplier;
end
end
