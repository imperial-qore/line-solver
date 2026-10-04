function res = gaver_stehfest_get_omega(n)
% RES = GAVER_STEHFEST_GET_OMEGA(N)
%
% Gaver-Stehfest weights. N is rounded DOWN to even. Twin of the native Python
% api.lti.gaver_stehfest_get_omega.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if mod(n,2) == 1
    n = n - 1;
end
res = zeros(1, n);
ndiv2 = n/2;
for k = 1:n
    val = ((-1.0)^(ndiv2 + k)) * log(2.0);
    sumVal = 0.0;
    for j = floor((k+1)/2):min(k, ndiv2)
        val2 = (j^(ndiv2 + 1));
        val2 = val2 / factorial(ndiv2);
        val2 = val2 * nchoosek(ndiv2, j);
        val2 = val2 * nchoosek(2*j, j);
        val2 = val2 * nchoosek(j, k-j);
        sumVal = sumVal + val2;
    end
    res(k) = val * sumVal;
end
end
