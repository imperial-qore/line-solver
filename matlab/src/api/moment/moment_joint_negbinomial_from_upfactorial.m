function bm = moment_joint_negbinomial_from_upfactorial(fp)
% bm = moment_joint_negbinomial_from_upfactorial(fp)
%
% Converts joint upward-factorial moments into the joint negative-binomial
% moments b-_(i_1,...,i_d) = E[prod_j nchoosek(N_j+i_j-1,i_j)] = f+_(i) /
% prod_j (i_j!).
%
% Input:
%   fp: array of size (n_1+1)x...x(n_d+1) holding the joint upward-factorial
%       moments, element (i_1+1,...,i_d+1) being the moment of
%       multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   bm: array of the same size holding the joint negative-binomial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

bm = moment_jointtrans(fp, 'negbinomial_from_upfactorial');
end
