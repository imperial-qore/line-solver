function m = moment_joint_raw_from_central(mc, mu)
% m = moment_joint_raw_from_central(mc, mu)
%
% Converts the joint central moments of a random vector into the joint power
% (raw) moments, by the multi-index binomial theorem,
%
%   m_(i) = sum_(k<=i) prod_j nchoosek(i_j,k_j) mu_j^(i_j-k_j) * mc_(k)
%
% The mean vector must be supplied separately, since the first-order central
% moments are zero and carry no information on it. Inverse of
% moment_joint_central_from_raw.
%
% Input:
%   mc: array of size (n_1+1)x...x(n_d+1) holding the joint central moments
%   mu: vector of length d holding the means E[N_1],...,E[N_d]
%
% Output:
%   m: array of the same size as mc holding the joint power moments
%
% Example:
%   m = moment_joint_raw_from_central(moment_joint_central_from_raw(m0), mu);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.

sz = moment_tensorsize(mc);
d = numel(sz);
if numel(mu) ~= d
    line_error(mfilename,'The mean vector mu must have one entry per dimension of mc.');
end
m = mc;
for mode = 1:d
    n = sz(mode)-1;
    T = zeros(n+1,n+1);
    for i = 0:n
        for k = 0:i
            T(i+1,k+1) = nchoosek(i,k) * mu(mode)^(i-k);
        end
    end
    m = moment_tensortrans(m, T, mode);
end
end
