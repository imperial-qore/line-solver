function f = moment_factorial_from_factcumulant(kappa)
% f = moment_factorial_from_factcumulant(kappa)
%
% Converts the factorial cumulants of a discrete random variable N into its
% factorial moments f_n = E[N(N-1)...(N-n+1)], by running the recursion
%
%   f_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * f_(n-k)
%
% forward with f_0 = 1. Inverse of moment_factcumulant_from_factorial.
%
% Input:
%   kappa: vector of length n+1 holding the factorial cumulants of order
%          0,...,n. Element 1 is ignored
%
% Output:
%   f: vector of length n+1 holding f_0,...,f_n, with the same orientation as
%      kappa and f(1) = 1
%
% Example:
%   f = moment_factorial_from_factcumulant([0, 2, 0, 0]);
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

f = moment_raw_from_cumulant(kappa);
end
