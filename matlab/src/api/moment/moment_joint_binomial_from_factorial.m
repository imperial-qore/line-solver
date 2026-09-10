function b = moment_joint_binomial_from_factorial(f)
% b = moment_joint_binomial_from_factorial(f)
%
% Converts joint factorial moments into the joint binomial moments
% b_(i_1,...,i_d) = E[prod_j nchoosek(N_j,i_j)] = f_(i) / prod_j (i_j!).
%
% Input:
%   f: array of size (n_1+1)x...x(n_d+1) holding the joint factorial
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   b: array of the same size holding the joint binomial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

b = moment_jointtrans(f, 'binomial_from_factorial');
end
