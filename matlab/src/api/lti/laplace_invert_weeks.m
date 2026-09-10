function f = laplace_invert_weeks(F, t, q, sigma, b)
% F = LAPLACE_INVERT_WEEKS(F, T, Q, SIGMA, B)
%
% Invert a Laplace transform by the Laguerre (Weeks) series
%
%     f(t) = sum_{n>=0} q_n l_n(t)
%
% of W. Weeks, "Numerical inversion of Laplace transforms using Laguerre
% functions", J. ACM 13, 1966, in the form used by P. G. Harrison and
% W. J. Knottenbelt, "Passage Time Distributions in Large Markov Chains",
% 2002, Sec. 4.1-4.3, after J. Abate, G. Choudhury and W. Whitt, "On the
% Laguerre method for numerically inverting Laplace transforms", INFORMS
% J. Computing 8(4), 1996.
%
% F     transform handle, called with a COMPLEX argument
% T     time points (vector); t <= 0 returns 0
% Q     Laguerre coefficients from LAPLACE_WEEKS_COEFFS; omit or leave empty
%       to have LAPLACE_WEEKS_SCALING pick sigma and b and compute them
% SIGMA, B  the exponential damping and scaling parameters that Q was built
%       with; required when Q is supplied
%
% Unlike Euler and Talbot the coefficients do not depend on t, so one
% coefficient set serves an arbitrary number of time points. That is the
% property this method is here for: the transform is evaluated 2*p0 times in
% total, not 2*p0 times per t.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(q)
    [sigma, b, q] = laplace_weeks_scaling(F);
elseif nargin < 5
    line_error(mfilename, 'When Q is supplied, SIGMA and B must be supplied with it: they are the parameters Q was built at, and evaluating Q at other values silently returns a different function.');
end

f = zeros(size(t));
nterms = laplace_weeks_nterms(q);
for i = 1:numel(t)
    if t(i) <= 0
        continue
    end
    l = laplace_weeks_functions(b * t(i), nterms);
    f(i) = exp(sigma * b * t(i)) * (q(1:nterms) * l);
end
end

function nterms = laplace_weeks_nterms(q)
% Truncate at the FIRST index where the coefficients have decayed, never the
% last. The quadrature divides by r^n with r = 0.1^(4/p0) < 1, so beyond the
% genuine decay the entries are rounding noise amplified by r^-n: at n = 2*p0
% that factor is 1e8, and a rule that scanned for the last entry above a
% threshold summed 1e-8 of pure noise (worst error on Exp(2) was 2.7e-09
% instead of 1.2e-13). p0 = numel(q)/2 is the hard cap, which is exactly the
% decay point LAPLACE_WEEKS_SCALING searched for.
p0 = floor(numel(q)/2);
nterms = p0;
for n = 2:(p0-1)
    if abs(q(n)) <= 1e-13 && abs(q(n+1)) <= 1e-13
        nterms = n;
        return
    end
end
end

function l = laplace_weeks_functions(t, N)
% l_n(t) = exp(-t/2) L_n(t) by the stable three-term recursion of Sec. 4.1.
l = zeros(N, 1);
l(1) = exp(-t/2);
if N > 1
    l(2) = (1 - t) * l(1);
end
for n = 2:(N-1)
    l(n+1) = ((2*n - 1 - t)/n) * l(n) - ((n - 1)/n) * l(n-1);
end
end
