function res = laplace_invert_euler(F, t, n)
% RES = LAPLACE_INVERT_EULER(F, T, N)
%
% Invert a Laplace transform at T by the Euler (Abate-Whitt) method. F is a
% handle called with a COMPLEX argument.
%
% N defaults to 41 and is rounded UP to odd. THE LARGER N IS NOT THE BETTER
% ONE: the weights carry a factor 10^((n-1)/6) against an ALTERNATING sum, so
% accuracy is a race between the series converging and the cancellation eating
% the mantissa. Worst relative error on F(s) = 2/(s+2) over
% t in {0.1, 0.5, 1, 2}:
%
%     n   =    11      21      31      41      51      71      99
%     err = 4.4e-3  2.1e-6  1.6e-9  1.6e-10 4.7e-8  1.7e-4  1.4e+0
%
% At 99 the result is 140 per cent wrong and negative at some t. In double
% precision 41 is the optimum; the native Python twin carries the same table.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(n)
    n = 41;
end
if mod(n,2) == 0
    n = n + 1;
end
if t <= 0
    line_error(mfilename, 'The Laplace inversion time point must be positive.');
end

alpha = euler_get_alpha(n);
omega = euler_get_omega(n);
res = 0.0;
for i = 1:n
    res = res + real(omega(i) * F(alpha(i)/t));
end
res = res / t;
end
