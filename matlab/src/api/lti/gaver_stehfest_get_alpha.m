function res = gaver_stehfest_get_alpha(n)
% RES = GAVER_STEHFEST_GET_ALPHA(N)
%
% Gaver-Stehfest nodes alpha_k = k log 2, k = 1..N. N is rounded DOWN to even
% (the method is defined for even N only); Euler rounds UP instead, and the
% two conventions must not be swapped. Twin of the native Python
% api.lti.gaver_stehfest_get_alpha.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if mod(n,2) == 1
    n = n - 1;
end
res = (1:n) * log(2.0);
end
