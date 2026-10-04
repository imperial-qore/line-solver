function m = moment_raw_from_upfactorial(fp)
% m = moment_raw_from_upfactorial(fp)
%
% Converts the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a
% discrete random variable N into the power (raw) moments m_n = E[N^n] by
% means of the signed Stirling numbers of the second kind,
%
%   m_n = sum_{k=0}^{n} (-1)^(n-k) * S(n,k) * f_k^+
%
% Input:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, i.e. fp(i) is the moment
%       of order i-1 and fp(1) = f_0^+ = 1
%
% Output:
%   m: vector of length n+1 holding m_0,...,m_n, with the same orientation
%      as fp
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% m = moment_raw_from_upfactorial(moment_upfactorial_from_raw([1,2,6,22]))

fpcol = fp(:);
n = length(fpcol)-1;
S = moment_stirling2(n);
T = zeros(n+1,n+1);
for i = 0:n
    for j = 0:i
        T(i+1,j+1) = (-1)^(i-j) * S(i+1,j+1);
    end
end
m = T * fpcol;
if isrow(fp)
    m = m.';
end
end
