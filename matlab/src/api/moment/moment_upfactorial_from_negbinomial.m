function fp = moment_upfactorial_from_negbinomial(bm)
% fp = moment_upfactorial_from_negbinomial(bm)
%
% Converts the negative-binomial moments b_n^- = E[nchoosek(N+n-1,n)] of a
% discrete random variable N into the upward-factorial moments
% f_n^+ = E[N(N+1)...(N+n-1)] via the one-to-one correspondence
%
%   f_n^+ = n! * b_n^-
%
% Input:
%   bm: vector of length n+1 holding b_0^-,...,b_n^-, i.e. bm(i) is the moment
%       of order i-1 and bm(1) = b_0^- = 1
%
% Output:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, with the same
%       orientation as bm
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (7).
%
% Example:
% fp = moment_upfactorial_from_negbinomial([1,2,4,22/3])

bmcol = bm(:);
n = length(bmcol)-1;
fp = zeros(n+1,1);
for i = 0:n
    fp(i+1) = factorial(i) * bmcol(i+1);
end
if isrow(bm)
    fp = fp.';
end
end
