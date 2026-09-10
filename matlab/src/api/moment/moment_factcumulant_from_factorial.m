function kappa = moment_factcumulant_from_factorial(f)
% kappa = moment_factcumulant_from_factorial(f)
%
% Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete random
% variable N into its factorial cumulants, the coefficients of the logarithm of
% the probability generating function expanded about z = 1,
%
%   log E[z^N] = sum_{n>=1} kappa_n (z-1)^n / n!
%
% The factorial cumulants stand to the factorial moments exactly as the
% cumulants stand to the power moments, so the same recursion applies,
%
%   f_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * f_(n-k)
%
% For a Poisson variable of rate lambda all factorial cumulants beyond the
% first vanish, which makes them the natural measure of departure from Poisson
% behaviour in the counting process of a MAP.
%
% Input:
%   f: vector of length n+1 holding f_0,...,f_n, i.e. f(i) is the moment of
%      order i-1 and f(1) = 1
%
% Output:
%   kappa: vector of length n+1 holding the factorial cumulants of order
%          0,...,n, with the same orientation as f and element 1 equal to 0
%
% Example:
%   kappa = moment_factcumulant_from_factorial([1, 2, 4, 8]);
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

kappa = moment_cumulant_from_raw(f);
end
