function q = laplace_weeks_coeffs(F, sigma, b, p0)
% Q = LAPLACE_WEEKS_COEFFS(F, SIGMA, B, P0)
%
% Laguerre coefficients q_n, n = 0..2*p0-1, of the damped and scaled function
%
%     f_{sigma,b}(t) = exp(-sigma t) f(t/b)
%
% whose Laguerre generating function is (Harrison and Knottenbelt 2002,
% Sec. 4.2)
%
%     Q_{sigma,b}(z) = b/(1-z) * L( b(1+z)/(2(1-z)) + b*sigma )
%
% and q_n = (1/(2 pi i)) contour-integral Q(z)/z^{n+1} dz on |z| = r.
%
% NOTE ON THE PAPER. Eq. 10 as printed carries the factor (1-z) rather than
% 1/(1-z). The scaled form quoted above, printed later in the same section,
% carries 1/(1-z) and is the correct one: with l_n(t) = exp(-t/2) L_n(t) the
% transform of l_n is (s-1/2)^n/(s+1/2)^{n+1}, so L(s) = Q(z)/(s+1/2) with
% z = (s-1/2)/(s+1/2), and s+1/2 = 1/(1-z). Implementing the printed (1-z)
% gives a wrong answer at every t (163 per cent at t = 0.1 on Exp(2)), so the
% discrepancy is a typo and not a convention.
%
% Sec. 4.3 fixes the number of trapezoids at 2*p0 and the radius at
% r = 0.1^(4/p0) for every n, rather than letting them grow with n. The
% resulting quadrature is a discrete Fourier transform of Q sampled on the
% circle, so all 2*p0 coefficients come out of one FFT and the transform is
% evaluated 2*p0 times in total.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(p0)
    p0 = 200;
end
if nargin < 3 || isempty(b)
    b = 1;
end
if nargin < 2 || isempty(sigma)
    sigma = 0;
end
if b <= 0
    line_error(mfilename, 'The scaling parameter b must be positive.');
end

N = 2 * p0;
r = 0.1^(4/p0);
j = (0:N-1)';
z = r * exp(2i * pi * j / N);
s = b * (1 + z) ./ (2 * (1 - z)) + b * sigma;

Qz = zeros(N, 1);
for k = 1:N
    Qz(k) = b / (1 - z(k)) * F(s(k));
end

q = real(fft(Qz)).' / N;
q = q ./ (r .^ (0:N-1));
end
