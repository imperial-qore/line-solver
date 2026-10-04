function kappa = moment_joint_factcumulant_from_factorial(f)
% kappa = moment_joint_factcumulant_from_factorial(f)
%
% Converts the joint factorial moments of a discrete random vector into its
% joint factorial cumulants, the coefficients of the logarithm of the joint
% probability generating function expanded about z = (1,...,1),
%
%   log E[prod_j z_j^(N_j)] = sum_(a ~= 0) kappa_a prod_j (z_j-1)^(a_j) / a_j!
%
% They stand to the joint factorial moments exactly as the joint cumulants
% stand to the joint power moments, so the same recursion applies. For a
% multivariate Poisson vector with independent components every joint factorial
% cumulant of order two or more vanishes; for the per-class counts of a marked
% MAP they measure the departure from independent Poisson marking.
%
% Input:
%   f: array of size (n_1+1)x...x(n_d+1) holding the joint factorial moments,
%      with element 1 equal to 1
%
% Output:
%   kappa: array of the same size holding the joint factorial cumulants,
%          element 1 being 0
%
% Example:
%   kappa = moment_joint_factcumulant_from_factorial(f);
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

kappa = moment_joint_cumulant_from_raw(f);
end
