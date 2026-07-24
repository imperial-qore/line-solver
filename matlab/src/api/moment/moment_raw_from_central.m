function m = moment_raw_from_central(mc, m1)
% m = moment_raw_from_central(mc, m1)
%
% Converts the central moments m_n^c = E[(N-m_1)^n] of a random variable N
% into the power (raw) moments m_n = E[N^n] by means of the inverse binomial
% transform in the variation that involves the mean m_1,
%
%   m_n = sum_{k=0}^{n} nchoosek(n,k) * m_k^c * m_1^(n-k)
%
% The mean must be supplied separately since m_1^c = 0 carries no information
% on it. The conversion also holds for continuous random variables.
%
% Input:
%   mc: vector of length n+1 holding m_0^c,...,m_n^c, i.e. mc(i) is the moment
%       of order i-1 and mc(1) = m_0^c = 1
%   m1: mean of N
%
% Output:
%   m: vector of length n+1 holding m_0,...,m_n, with the same orientation
%      as mc
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% m = moment_raw_from_central(moment_central_from_raw([1,2,6,22]), 2)

mccol = mc(:);
n = length(mccol)-1;
if ~isscalar(m1)
    line_error(mfilename,'The mean m1 must be a scalar.');
end
m = zeros(n+1,1);
for i = 0:n
    for k = 0:i
        m(i+1) = m(i+1) + nchoosek(i,k) * mccol(k+1) * m1^(i-k);
    end
end
if isrow(mc)
    m = m.';
end
end
