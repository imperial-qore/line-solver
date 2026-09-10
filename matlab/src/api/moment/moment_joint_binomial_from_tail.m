function b = moment_joint_binomial_from_tail(t)
% b = moment_joint_binomial_from_tail(t)
%
% Converts the joint survival probabilities t_(m_1,...,m_d) = P(N_1 >= m_1,
% ..., N_d >= m_d) of a nonnegative integer random vector into its joint
% binomial moments,
%
%   b_(k) = E[prod_j nchoosek(N_j,k_j)]
%         = sum_(m>=k) prod_j nchoosek(m_j-1,k_j-1) * t_(m)
%
% The transform is the tensor product of the univariate one, which is what
% makes the mode-by-mode evaluation legitimate: an entry with k_j = 0 selects
% m_j = 0, and t_(0,m_2,...) is by construction the marginal survival array of
% the remaining coordinates. As in the univariate case it is upper triangular,
% so the array must cover the joint support to be exact; truncating gives lower
% bounds.
%
% Input:
%   t: array of size (n_1+1)x...x(n_d+1) holding the joint survival
%      probabilities, element 1 being 1
%
% Output:
%   b: array of the same size holding the joint binomial moments
%
% Example:
%   b = moment_joint_binomial_from_tail(t);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

b = moment_jointtrans(t, 'binomial_from_tail');
end
