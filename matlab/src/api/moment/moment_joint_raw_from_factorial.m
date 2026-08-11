function m = moment_joint_raw_from_factorial(f)
% m = moment_joint_raw_from_factorial(f)
%
% Converts joint factorial moments into joint power (raw) moments, by the
% Stirling numbers of the second kind along every dimension. Inverse of
% moment_joint_factorial_from_raw.
%
% Input:
%   f: array of size (n_1+1)x...x(n_d+1) holding the joint factorial
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   m: array of the same size holding the joint raw moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

m = moment_jointtrans(f, 'raw_from_factorial');
end
