function x = moment_binotransinv(y)
% x = moment_binotransinv(y)
%
% Inverse binomial transform of the sequence y_0,y_1,...,y_n,
%
%   x_n = sum_{k=0}^{n} nchoosek(n,k) * y_k
%
% Input:
%   y: vector of length n+1 holding y_0,...,y_n, i.e. y(i) is the element of
%      order i-1
%
% Output:
%   x: vector of length n+1 holding x_0,...,x_n, with the same orientation
%      as y
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (9).
%
% Example:
% x = moment_binotransinv(moment_binotrans([1,2,5,15]))

ycol = y(:);
n = length(ycol)-1;
x = zeros(n+1,1);
for i = 0:n
    for k = 0:i
        x(i+1) = x(i+1) + nchoosek(i,k) * ycol(k+1);
    end
end
if isrow(y)
    x = x.';
end
end
