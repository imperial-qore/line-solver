function [sigma, b, q] = laplace_weeks_scaling(F, p0, tol)
% [SIGMA, B, Q] = LAPLACE_WEEKS_SCALING(F, P0, TOL)
%
% Automatic selection of the exponential damping SIGMA and the scaling B for
% the Laguerre inversion, following the algorithm of Fig. 1 of P. G. Harrison
% and W. J. Knottenbelt, "Passage Time Distributions in Large Markov Chains",
% 2002. The search accepts the first (sigma, b) at which the coefficients have
% decayed by term P0:
%
%     sigma = 0; b = 1
%     while |q_p0| > tol or |q_{p0+1}| > tol
%         sigma = 0.001 if sigma == 0 else 2*sigma
%         if sigma > 0.2
%             b = b + 4;  if b > 10, no suitable parameters exist
%             sigma = 0
%
% Raising b too far is numerically counterproductive and excessive damping is
% unstable in finite precision, which is why the search is bounded rather than
% unbounded. When it exhausts the box the failure is REFUSED BY NAME: a
% density with a discontinuity in itself or its derivatives has no usable
% Laguerre representation (Sec. 4.2), and returning the last iterate would
% report noise as an answer. Use 'euler' for those, at roughly 50 transform
% evaluations per time point.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(p0)
    p0 = 200;
end
if nargin < 3 || isempty(tol)
    tol = 1e-10;
end

sigma = 0;
b = 1;
while true
    q = laplace_weeks_coeffs(F, sigma, b, p0);
    if abs(q(p0+1)) <= tol && abs(q(p0+2)) <= tol
        return
    end
    if sigma == 0
        sigma = 0.001;
    else
        sigma = 2 * sigma;
    end
    if sigma > 0.2
        b = b + 4;
        if b > 10
            line_error(mfilename, 'no suitable scaling parameters were found for the Laguerre inversion: the transform''s density is not smooth enough for a Laguerre series. Use laplace_invert(F,t,''euler'') instead.');
        end
        sigma = 0;
    end
end
end
