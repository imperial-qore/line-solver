function res = euler_get_eta(n)
% RES = EULER_GET_ETA(N)
%
% Euler (Abate-Whitt) weights eta_i, i = 1..N, N odd. Twin of the native
% Python api.lti.euler_get_eta.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

res = zeros(1, n);
res(1) = 0.5;
for i = 2:((n+1)/2)
    res(i) = 1.0;
end
res(n) = 1.0 / (2.0^((n-1)/2.0));
for i = 1:(((n-1)/2) - 1)
    res(n-i) = res(n-i+1) + (2.0^((1-n)/2.0)) * nchoosek((n-1)/2, i);
end
end
