function arr = talbot_get_alpha(n)
% ARR = TALBOT_GET_ALPHA(N)
%
% Talbot contour nodes. Twin of the native Python api.lti.talbot_get_alpha.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

arr = zeros(1, n);
arr(1) = complex(2.0 * n / 5.0, 0.0);
for i = 2:n
    theta = (i-1) * pi / n;
    arr(i) = complex(2*(i-1)*pi/5 * (1.0/tan(theta)), 2*(i-1)*pi/5);
end
end
