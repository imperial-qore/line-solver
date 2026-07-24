function bm = moment_negbinomial_from_binomial(b)
% bm = moment_negbinomial_from_binomial(b)
%
% Converts the binomial moments b_n = E[nchoosek(N,n)] of a discrete random
% variable N into the negative-binomial moments b_n^- = E[nchoosek(N+n-1,n)]
% by means of the shifted binomial transform
%
%   b_n^- = sum_{k=1}^{n} nchoosek(n-1,k-1) * b_k   for n >= 1
%   b_0^- = 1
%
% Input:
%   b: vector of length n+1 holding b_0,...,b_n, i.e. b(i) is the moment of
%      order i-1 and b(1) = b_0 = 1
%
% Output:
%   bm: vector of length n+1 holding b_0^-,...,b_n^-, with the same
%       orientation as b
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (14).
%
% Example:
% bm = moment_negbinomial_from_binomial([1,2,2,4/3])

bcol = b(:);
n = length(bcol)-1;
bm = zeros(n+1,1);
bm(1) = 1;
for i = 1:n
    for k = 1:i
        bm(i+1) = bm(i+1) + nchoosek(i-1,k-1) * bcol(k+1);
    end
end
if isrow(b)
    bm = bm.';
end
end
