function bm = moment_negbinomial_from_upfactorial(fp)
% bm = moment_negbinomial_from_upfactorial(fp)
%
% Converts the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a
% discrete random variable N into the negative-binomial moments
% b_n^- = E[nchoosek(N+n-1,n)] via the one-to-one correspondence
%
%   b_n^- = f_n^+ / n!
%
% Input:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, i.e. fp(i) is the moment
%       of order i-1 and fp(1) = f_0^+ = 1
%
% Output:
%   bm: vector of length n+1 holding b_0^-,...,b_n^-, with the same
%       orientation as fp
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, eq. (7).
%
% Example:
% bm = moment_negbinomial_from_upfactorial([1,2,8,44])

fpcol = fp(:);
n = length(fpcol)-1;
bm = zeros(n+1,1);
for i = 0:n
    bm(i+1) = fpcol(i+1) / factorial(i);
end
if isrow(fp)
    bm = bm.';
end
end
