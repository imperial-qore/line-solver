function mc = moment_central_from_raw(m)
% mc = moment_central_from_raw(m)
%
% Converts the power (raw) moments m_n = E[N^n] of a random variable N into
% the central moments m_n^c = E[(N-m_1)^n] by means of the binomial transform
% in the variation that involves the mean m_1,
%
%   m_n^c = sum_{k=0}^{n} (-1)^(n-k) * nchoosek(n,k) * m_k * m_1^(n-k)
%
% The conversion also holds for continuous random variables.
%
% Input:
%   m: vector of length n+1 holding m_0,...,m_n, i.e. m(i) is the moment of
%      order i-1 and m(1) = m_0 = 1. At least the mean m_1 must be given,
%      hence n >= 1
%
% Output:
%   mc: vector of length n+1 holding m_0^c,...,m_n^c, with the same
%       orientation as m. By construction m_0^c = 1 and m_1^c = 0
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% mc = moment_central_from_raw([1,2,6,22])

mcol = m(:);
n = length(mcol)-1;
if n < 1
    line_error(mfilename,'The mean m_1 is required for this conversion, hence m must have at least 2 elements.');
end
m1 = mcol(2);
mc = zeros(n+1,1);
for i = 0:n
    for k = 0:i
        mc(i+1) = mc(i+1) + (-1)^(i-k) * nchoosek(i,k) * mcol(k+1) * m1^(i-k);
    end
end
if isrow(m)
    mc = mc.';
end
end
