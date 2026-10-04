function fp = moment_joint_upfactorial_from_raw(m)
% fp = moment_joint_upfactorial_from_raw(m)
%
% Converts joint power (raw) moments into the joint upward-factorial moments
% f+_(i_1,...,i_d) = E[prod_j N_j(N_j+1)...(N_j+i_j-1)], by the Stirling
% cycle numbers along every dimension.
%
% Input:
%   m: array of size (n_1+1)x...x(n_d+1) holding the joint raw
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   fp: array of the same size holding the joint upward-factorial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

fp = moment_jointtrans(m, 'upfactorial_from_raw');
end
