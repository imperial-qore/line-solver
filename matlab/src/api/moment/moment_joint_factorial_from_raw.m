function f = moment_joint_factorial_from_raw(m)
% f = moment_joint_factorial_from_raw(m)
%
% Converts the joint power (raw) moments m_(i_1,...,i_d) = E[prod_j N_j^(i_j)]
% of a random vector (N_1,...,N_d) into the joint factorial moments
% f_(i_1,...,i_d) = E[prod_j (N_j)_(i_j)], where (N)_i = N(N-1)...(N-i+1), by
% applying the signed Stirling numbers of the first kind separately along
% every dimension,
%
%   f_(i) = sum_(k) prod_j s(i_j,k_j) * m_(k)
%
% The joint conversion is the Kronecker product of the univariate ones,
% which is what makes the mode-by-mode evaluation legitimate. Only the
% cumulant and the central conversions are not of this separable form.
%
% Input:
%   m: array of size (n_1+1)x...x(n_d+1) holding the joint raw
%      moments, element (i_1+1,...,i_d+1) being the moment of
%      multi-order (i_1,...,i_d) and element 1 being 1
%
% Output:
%   f: array of the same size holding the joint factorial moments
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

f = moment_jointtrans(m, 'factorial_from_raw');
end
