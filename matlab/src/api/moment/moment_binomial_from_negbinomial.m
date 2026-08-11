function b = moment_binomial_from_negbinomial(bm)
% b = moment_binomial_from_negbinomial(bm)
%
% Converts the negative-binomial moments b_n^- = E[nchoosek(N+n-1,n)] of a
% discrete random variable N into the binomial moments b_n = E[nchoosek(N,n)]
% by means of the shifted binomial transform
%
%   b_n = sum_{k=1}^{n} (-1)^(n-k) * nchoosek(n-1,k-1) * b_k^-   for n >= 1
%   b_0 = 1
%
% Input:
%   bm: vector of length n+1 holding b_0^-,...,b_n^-, i.e. bm(i) is the moment
%       of order i-1 and bm(1) = b_0^- = 1
%
% Output:
%   b: vector of length n+1 holding b_0,...,b_n, with the same orientation
%      as bm
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (14).
%
% Example:
% b = moment_binomial_from_negbinomial([1,2,4,22/3])

bmcol = bm(:);
n = length(bmcol)-1;
b = zeros(n+1,1);
b(1) = 1;
for i = 1:n
    for k = 1:i
        b(i+1) = b(i+1) + (-1)^(i-k) * nchoosek(i-1,k-1) * bmcol(k+1);
    end
end
if isrow(bm)
    b = b.';
end
end
