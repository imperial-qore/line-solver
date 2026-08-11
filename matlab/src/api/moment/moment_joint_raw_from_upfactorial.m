function m = moment_joint_raw_from_upfactorial(fp)
% m = moment_joint_raw_from_upfactorial(fp)
%
% Converts joint upward-factorial moments into joint power (raw) moments, by
% the signed Stirling numbers of the second kind along every dimension.
% Inverse of moment_joint_upfactorial_from_raw.
%
% Input:
%   fp: array of size (n_1+1)x...x(n_d+1) holding the joint upward-factorial
%       moments, element (i_1+1,...,i_d+1) being the moment of
%       multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   m: array of the same size holding the joint raw moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

m = moment_jointtrans(fp, 'raw_from_upfactorial');
end
