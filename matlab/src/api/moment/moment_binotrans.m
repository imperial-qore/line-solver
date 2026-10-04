function y = moment_binotrans(x)
% y = moment_binotrans(x)
%
% Binomial transform of the sequence x_0,x_1,...,x_n into y_0,y_1,...,y_n,
%
%   y_n = sum_{k=0}^{n} (-1)^(n-k) * nchoosek(n,k) * x_k
%
% Applied to a moment sequence m_i = E[X^i] it returns the moments of the
% unit downshift, y_i = E[(X-1)^i]. It is not an involution: its inverse is
% moment_binotransinv, the unsigned transform.
%
% Input:
%   x: vector of length n+1 holding x_0,...,x_n, i.e. x(i) is the element of
%      order i-1
%
% Output:
%   y: vector of length n+1 holding y_0,...,y_n, with the same orientation
%      as x
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (8).
%
% Example:
% y = moment_binotrans([1,2,5,15])

xcol = x(:);
n = length(xcol)-1;
y = zeros(n+1,1);
for i = 0:n
    for k = 0:i
        y(i+1) = y(i+1) + (-1)^(i-k) * nchoosek(i,k) * xcol(k+1);
    end
end
if isrow(x)
    y = y.';
end
end
