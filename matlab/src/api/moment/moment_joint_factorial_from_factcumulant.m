function f = moment_joint_factorial_from_factcumulant(kappa)
% f = moment_joint_factorial_from_factcumulant(kappa)
%
% Converts the joint factorial cumulants of a discrete random vector into its
% joint factorial moments. Inverse of
% moment_joint_factcumulant_from_factorial.
%
% Input:
%   kappa: array of size (n_1+1)x...x(n_d+1) holding the joint factorial
%          cumulants; element 1 is ignored
%
% Output:
%   f: array of the same size holding the joint factorial moments, element 1
%      being 1
%
% Example:
%   f = moment_joint_factorial_from_factcumulant(kappa);
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

f = moment_joint_raw_from_cumulant(kappa);
end
