function sigma = moment_stirlingcycle(n)
% sigma = moment_stirlingcycle(n)
%
% Triangle of the Stirling cycle numbers (unsigned Stirling numbers of the
% first kind), sigma(i,j) = (-1)^(i-j) * s(i,j), obtained from the recursion
%
%   sigma(i,j) = (i-1)*sigma(i-1,j) + sigma(i-1,j-1)   for j > 0
%   sigma(0,0) = 1,   sigma(i,0) = 0 for i > 0
%
% These numbers are the coefficients that convert power moments into
% upward-factorial moments.
%
% Input:
%   n: maximum order (n >= 0)
%
% Output:
%   sigma: (n+1)x(n+1) matrix with sigma(i+1,j+1) = sigma(i,j) in the 0-based
%          notation of the reference. Entries with j > i are zero.
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (12).
%
% Example:
% sigma = moment_stirlingcycle(3)

if ~isscalar(n) || n < 0 || n ~= round(n)
    line_error(mfilename,'The maximum order n must be a nonnegative integer.');
end
sigma = zeros(n+1,n+1);
sigma(1,1) = 1;
for i = 1:n
    for j = 1:i
        sigma(i+1,j+1) = (i-1)*sigma(i,j+1) + sigma(i,j);
    end
end
end
