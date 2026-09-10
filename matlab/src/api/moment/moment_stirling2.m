function S = moment_stirling2(n)
% S = moment_stirling2(n)
%
% Triangle of the Stirling numbers of the second kind S(i,j), implicitly
% defined by the expansion of a power into falling factorials
%
%   x^i = sum_{j=0}^{i} S(i,j) x(x-1)(x-2)...(x-j+1)
%
% and computed from the recursion
%
%   S(i,j) = j*S(i-1,j) + S(i-1,j-1)   for j > 0
%   S(0,0) = 1,   S(i,0) = 0 for i > 0
%
% These numbers are the coefficients that convert factorial moments back into
% power moments.
%
% Input:
%   n: maximum order (n >= 0)
%
% Output:
%   S: (n+1)x(n+1) matrix with S(i+1,j+1) = S(i,j) in the 0-based notation of
%      the reference. Entries with j > i are zero.
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (11).
%
% Example:
% S = moment_stirling2(3)

if ~isscalar(n) || n < 0 || n ~= round(n)
    line_error(mfilename,'The maximum order n must be a nonnegative integer.');
end
S = zeros(n+1,n+1);
S(1,1) = 1;
for i = 1:n
    for j = 1:i
        S(i+1,j+1) = j*S(i,j+1) + S(i,j);
    end
end
end
