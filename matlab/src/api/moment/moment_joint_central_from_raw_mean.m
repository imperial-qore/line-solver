function mc = moment_joint_central_from_raw_mean(m, mu)
% mc = moment_joint_central_from_raw_mean(m, mu)
%
% Converts the joint power (raw) moments of a random vector into the joint
% central moments about a given mean vector. Same conversion as
% moment_joint_central_from_raw, with the means supplied rather than read off
% the array, so that it also applies when the array does not carry the
% first-order entries.
%
% Input:
%   m: array of size (n_1+1)x...x(n_d+1) holding the joint power moments
%   mu: vector of length d holding the means E[N_1],...,E[N_d]
%
% Output:
%   mc: array of the same size as m holding the joint central moments
%
% Example:
%   mc = moment_joint_central_from_raw_mean(m, [1.5, 2.5]);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.

sz = moment_tensorsize(m);
d = numel(sz);
if numel(mu) ~= d
    line_error(mfilename,'The mean vector mu must have one entry per dimension of m.');
end
mc = m;
for mode = 1:d
    n = sz(mode)-1;
    T = zeros(n+1,n+1);
    for i = 0:n
        for k = 0:i
            T(i+1,k+1) = nchoosek(i,k) * (-mu(mode))^(i-k);
        end
    end
    mc = moment_tensortrans(mc, T, mode);
end
end
