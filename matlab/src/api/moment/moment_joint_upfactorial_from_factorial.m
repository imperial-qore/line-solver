function fp = moment_joint_upfactorial_from_factorial(f)
% fp = moment_joint_upfactorial_from_factorial(f)
%
% Converts joint factorial moments into joint upward-factorial moments, by
% the Lah numbers along every dimension.
%
% Input:
%   f: array of size (n_1+1)x...x(n_d+1) holding the joint factorial
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   fp: array of the same size holding the joint upward-factorial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

fp = moment_jointtrans(f, 'upfactorial_from_factorial');
end
