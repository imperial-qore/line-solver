function L = moment_lah(n)
% L = moment_lah(n)
%
% Triangle of the Lah numbers L(i,j) = (i!/j!)*nchoosek(i-1,j-1), which link
% the factorial moments to the upward-factorial moments. The triangle is built
% from the equivalent recursion
%
%   L(i,j) = L(i-1,j-1) + (i+j-1)*L(i-1,j)   for j > 0
%   L(0,0) = 1,   L(i,0) = 0 for i > 0
%
% which avoids the overflow of the explicit factorial form for large orders.
%
% Input:
%   n: maximum order (n >= 0)
%
% Output:
%   L: (n+1)x(n+1) matrix with L(i+1,j+1) = L(i,j) in the 0-based notation of
%      the reference. Entries with j > i are zero.
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
% I. Lah. Eine neue Art von Zahlen, ihre Eigenschaften und Anwendung in der
% mathematischen Statistik. Mitteilungsbl. Math. Statist., 7:203-212, 1955.
%
% Example:
% L = moment_lah(3)

if ~isscalar(n) || n < 0 || n ~= round(n)
    line_error(mfilename,'The maximum order n must be a nonnegative integer.');
end
L = zeros(n+1,n+1);
L(1,1) = 1;
for i = 1:n
    for j = 1:i
        L(i+1,j+1) = L(i,j) + (i+j-1)*L(i,j+1);
    end
end
end
