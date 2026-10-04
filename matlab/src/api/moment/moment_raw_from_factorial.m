function m = moment_raw_from_factorial(f)
% m = moment_raw_from_factorial(f)
%
% Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete
% random variable N into the power (raw) moments m_n = E[N^n] by means of the
% Stirling numbers of the second kind,
%
%   m_n = sum_{k=0}^{n} S(n,k) * f_k
%
% Input:
%   f: vector of length n+1 holding f_0,...,f_n, i.e. f(i) is the moment of
%      order i-1 and f(1) = f_0 = 1
%
% Output:
%   m: vector of length n+1 holding m_0,...,m_n, with the same orientation
%      as f
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (13).
%
% Example:
% m = moment_raw_from_factorial(moment_factorial_from_raw([1,2,6,22]))

fcol = f(:);
n = length(fcol)-1;
m = moment_stirling2(n) * fcol;
if isrow(f)
    m = m.';
end
end
