function f = moment_factorial_from_raw(m)
% f = moment_factorial_from_raw(m)
%
% Converts the power (raw) moments m_n = E[N^n] of a discrete random variable
% N into the factorial moments f_n = E[N(N-1)...(N-n+1)] by means of the
% signed Stirling numbers of the first kind,
%
%   f_n = sum_{k=0}^{n} s(n,k) * m_k
%
% Input:
%   m: vector of length n+1 holding m_0,...,m_n, i.e. m(i) is the moment of
%      order i-1 and m(1) = m_0 = 1
%
% Output:
%   f: vector of length n+1 holding f_0,...,f_n, with the same orientation
%      as m
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (13).
%
% Example:
% f = moment_factorial_from_raw([1,2,6,22])

mcol = m(:);
n = length(mcol)-1;
f = moment_stirling1(n) * mcol;
if isrow(m)
    f = f.';
end
end
