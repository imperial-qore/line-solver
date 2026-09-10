function f = moment_joint_factorial_from_upfactorial(fp)
% f = moment_joint_factorial_from_upfactorial(fp)
%
% Converts joint upward-factorial moments into joint factorial moments, by
% the signed Lah numbers along every dimension.
%
% Input:
%   fp: array of size (n_1+1)x...x(n_d+1) holding the joint upward-factorial
%       moments, element (i_1+1,...,i_d+1) being the moment of
%       multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   f: array of the same size holding the joint factorial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

f = moment_jointtrans(fp, 'factorial_from_upfactorial');
end
