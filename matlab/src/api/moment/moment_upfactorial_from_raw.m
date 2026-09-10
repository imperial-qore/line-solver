function fp = moment_upfactorial_from_raw(m)
% fp = moment_upfactorial_from_raw(m)
%
% Converts the power (raw) moments m_n = E[N^n] of a discrete random variable
% N into the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] by means of
% the Stirling cycle numbers,
%
%   f_n^+ = sum_{k=0}^{n} sigma(n,k) * m_k
%
% Upward-factorial moments are of use in moment-matching techniques for
% matrix-geometric and discrete phase-type distributions.
%
% Input:
%   m: vector of length n+1 holding m_0,...,m_n, i.e. m(i) is the moment of
%      order i-1 and m(1) = m_0 = 1
%
% Output:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, with the same
%       orientation as m
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% fp = moment_upfactorial_from_raw([1,2,6,22])

mcol = m(:);
n = length(mcol)-1;
fp = moment_stirlingcycle(n) * mcol;
if isrow(m)
    fp = fp.';
end
end
