function b = moment_joint_binomial_from_negbinomial(bm)
% b = moment_joint_binomial_from_negbinomial(bm)
%
% Converts joint negative-binomial moments into joint binomial moments, by
% the signed shifted binomial transform along every dimension. Inverse of
% moment_joint_negbinomial_from_binomial.
%
% Input:
%   bm: array of size (n_1+1)x...x(n_d+1) holding the joint negative-binomial
%       moments, element (i_1+1,...,i_d+1) being the moment of
%       multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   b: array of the same size holding the joint binomial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

b = moment_jointtrans(bm, 'binomial_from_negbinomial');
end
