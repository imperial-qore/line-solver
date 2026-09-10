function mc = moment_joint_central_from_tail(t)
% mc = moment_joint_central_from_tail(t)
%
% Joint central moments of a nonnegative integer random vector from its joint
% survival array, composing the four edges that separate the two vertices:
% tail -> binomial -> factorial -> raw -> central, with the means read off the
% raw array.
%
% This is the whole path from a solver that produces survival probabilities (a
% closed queueing network through its normalizing constants, a CTMC through its
% stationary distribution, a simulator through a histogram) to the covariances
% and the higher central moments.
%
% Input:
%   t: array of size (n_1+1)x...x(n_d+1) holding the joint survival
%      probabilities and covering the support, element 1 being 1
%
% Output:
%   mc: array of the same size holding the joint central moments. The entry of
%       multi-order e_j+e_l is the covariance of N_j and N_l
%
% Example:
%   mc = moment_joint_central_from_tail(t);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.

b = moment_joint_binomial_from_tail(t);
f = moment_joint_factorial_from_binomial(b);
mc = moment_joint_central_from_raw(moment_joint_raw_from_factorial(f));
end
