function mc = moment_joint_central_from_raw(m)
% mc = moment_joint_central_from_raw(m)
%
% Converts the joint power (raw) moments of a random vector (N_1,...,N_d) into
% the joint central moments mc_(i_1,...,i_d) = E[prod_j (N_j - E N_j)^(i_j)].
%
% The conversion is the multi-index binomial theorem, which is again separable
% but with a different shift per dimension,
%
%   mc_(i) = sum_(k<=i) prod_j (-1)^(i_j-k_j) nchoosek(i_j,k_j) mu_j^(i_j-k_j)
%            * m_(k)
%
% The means mu_j = m_(e_j) are read off the array itself, so every dimension
% must carry at least the first order. The entry of multi-order e_j+e_l is the
% covariance of N_j and N_l.
%
% Input:
%   m: array of size (n_1+1)x...x(n_d+1) holding the joint power moments, with
%      every n_j >= 1
%
% Output:
%   mc: array of the same size holding the joint central moments
%
% Example:
%   mc = moment_joint_central_from_raw(m);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.

sz = moment_tensorsize(m);
d = numel(sz);
if any(sz < 2)
    line_error(mfilename,'The means m_(e_j) are required for this conversion, hence every dimension of m must have at least 2 elements.');
end
mu = zeros(1,d);
stride = cumprod([1, sz(1:end-1)]);
mv = reshape(m, [], 1);
for j = 1:d
    mu(j) = mv(1 + stride(j));
end
mc = moment_joint_central_from_raw_mean(m, mu);
end
