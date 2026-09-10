function kappa = moment_cumulant_from_raw(m)
% kappa = moment_cumulant_from_raw(m)
%
% Converts the power (raw) moments m_n = E[X^n] of a random variable X into its
% cumulants kappa_n, the coefficients of the cumulant generating function
% log E[exp(sX)] = sum_{n>=1} kappa_n s^n / n!.
%
% The conversion inverts the exponential-formula recursion
%
%   m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)
%
% equivalently kappa_n = sum_{pi in P(n)} (|pi|-1)! (-1)^(|pi|-1) prod_{B in pi}
% m_|B| over the set partitions of {1,...,n}. The first cumulants are
% kappa_1 = m_1, kappa_2 = m_2 - m_1^2 (the variance) and kappa_3 = m_3 -
% 3 m_1 m_2 + 2 m_1^3 (the third central moment). The conversion is not
% restricted to discrete random variables.
%
% Input:
%   m: vector of length n+1 holding m_0,...,m_n, i.e. m(i) is the moment of
%      order i-1 and m(1) = 1
%
% Output:
%   kappa: vector of length n+1 holding kappa_0,...,kappa_n, with the same
%          orientation as m. Element 1 is kappa_0 = 0, the value of the
%          cumulant generating function at the origin, and not m_0 = 1
%
% Example:
%   kappa = moment_cumulant_from_raw([1, 2, 6, 22]);
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

mcol = m(:);
n = length(mcol)-1;
kappa = zeros(n+1,1);
for i = 1:n
    acc = 0;
    for k = 1:(i-1)
        acc = acc + nchoosek(i-1,k-1) * kappa(k+1) * mcol(i-k+1);
    end
    kappa(i+1) = mcol(i+1) - acc;
end
if isrow(m)
    kappa = kappa.';
end
end
