function s = moment_stirling1(n)
% s = moment_stirling1(n)
%
% Triangle of the signed Stirling numbers of the first kind s(i,j), defined as
% the coefficients of x^j in the falling factorial
%
%   sum_{j=0}^{i} s(i,j) x^j = x(x-1)(x-2)...(x-i+1)
%
% These numbers are the coefficients that convert power moments into factorial
% moments. They relate to the Stirling cycle numbers via
% s(i,j) = (-1)^(i-j) * sigma(i,j).
%
% Input:
%   n: maximum order (n >= 0)
%
% Output:
%   s: (n+1)x(n+1) matrix with s(i+1,j+1) = s(i,j) in the 0-based notation of
%      the reference. Entries with j > i are zero.
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (10) and eq. (12).
%
% Example:
% s = moment_stirling1(3)

sigma = moment_stirlingcycle(n);
s = zeros(n+1,n+1);
for i = 0:n
    for j = 0:i
        s(i+1,j+1) = (-1)^(i-j) * sigma(i+1,j+1);
    end
end
end
