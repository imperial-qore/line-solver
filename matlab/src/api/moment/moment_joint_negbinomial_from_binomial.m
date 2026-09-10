function bm = moment_joint_negbinomial_from_binomial(b)
% bm = moment_joint_negbinomial_from_binomial(b)
%
% Converts joint binomial moments into joint negative-binomial moments, by
% the shifted binomial transform nchoosek(i-1,k-1) along every dimension.
%
% Input:
%   b: array of size (n_1+1)x...x(n_d+1) holding the joint binomial
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   bm: array of the same size holding the joint negative-binomial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

bm = moment_jointtrans(b, 'negbinomial_from_binomial');
end
