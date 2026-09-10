function f = moment_factorial_from_upfactorial(fp)
% f = moment_factorial_from_upfactorial(fp)
%
% Converts the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a
% discrete random variable N into the factorial moments
% f_n = E[N(N-1)...(N-n+1)] by means of the Lah numbers,
%
%   f_n = sum_{k=1}^{n} (-1)^(n-k) * L(n,k) * f_k^+   for n >= 1
%   f_0 = 1
%
% Input:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, i.e. fp(i) is the moment
%       of order i-1 and fp(1) = f_0^+ = 1
%
% Output:
%   f: vector of length n+1 holding f_0,...,f_n, with the same orientation
%      as fp
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% f = moment_factorial_from_upfactorial([1,2,8,44])

fpcol = fp(:);
n = length(fpcol)-1;
L = moment_lah(n);
f = zeros(n+1,1);
f(1) = 1;
for i = 1:n
    for k = 1:i
        f(i+1) = f(i+1) + (-1)^(i-k) * L(i+1,k+1) * fpcol(k+1);
    end
end
if isrow(fp)
    f = f.';
end
end
