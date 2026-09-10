function m = moment_raw_from_cumulant(kappa)
% m = moment_raw_from_cumulant(kappa)
%
% Converts the cumulants kappa_n of a random variable X into its power (raw)
% moments m_n = E[X^n], by running the exponential-formula recursion forward,
%
%   m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)
%
% with m_0 = 1. Equivalently m_n = sum_{pi in P(n)} prod_{B in pi} kappa_|B|
% over the set partitions of {1,...,n}. Inverse of moment_cumulant_from_raw.
%
% Input:
%   kappa: vector of length n+1 holding kappa_0,...,kappa_n, i.e. kappa(i) is
%          the cumulant of order i-1. Element 1 is ignored, since kappa_0 = 0
%          carries no information
%
% Output:
%   m: vector of length n+1 holding m_0,...,m_n, with the same orientation as
%      kappa and m(1) = 1
%
% Example:
%   m = moment_raw_from_cumulant(moment_cumulant_from_raw([1, 2, 6, 22]));
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

kcol = kappa(:);
n = length(kcol)-1;
m = zeros(n+1,1);
m(1) = 1;
for i = 1:n
    acc = 0;
    for k = 1:i
        acc = acc + nchoosek(i-1,k-1) * kcol(k+1) * m(i-k+1);
    end
    m(i+1) = acc;
end
if isrow(kappa)
    m = m.';
end
end
