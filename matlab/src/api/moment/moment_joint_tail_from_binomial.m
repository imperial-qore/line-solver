function t = moment_joint_tail_from_binomial(b)
% t = moment_joint_tail_from_binomial(b)
%
% Converts the joint binomial moments of a nonnegative integer random vector
% into its joint survival probabilities. Inverse of
% moment_joint_binomial_from_tail.
%
% Input:
%   b: array of size (n_1+1)x...x(n_d+1) holding the joint binomial moments
%
% Output:
%   t: array of the same size holding the joint survival probabilities
%
% Example:
%   t = moment_joint_tail_from_binomial(moment_joint_binomial_from_tail(t0));
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

t = moment_jointtrans(b, 'tail_from_binomial');
end
