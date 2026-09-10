function fp = moment_upfactorial_from_factorial(f)
% fp = moment_upfactorial_from_factorial(f)
%
% Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete
% random variable N into the upward-factorial moments
% f_n^+ = E[N(N+1)...(N+n-1)] by means of the Lah numbers,
%
%   f_n^+ = sum_{k=1}^{n} L(n,k) * f_k   for n >= 1
%   f_0^+ = 1
%
% Input:
%   f: vector of length n+1 holding f_0,...,f_n, i.e. f(i) is the moment of
%      order i-1 and f(1) = f_0 = 1
%
% Output:
%   fp: vector of length n+1 holding f_0^+,...,f_n^+, with the same
%       orientation as f
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% fp = moment_upfactorial_from_factorial([1,2,4,8])

fcol = f(:);
n = length(fcol)-1;
L = moment_lah(n);
fp = zeros(n+1,1);
fp(1) = 1;
for i = 1:n
    for k = 1:i
        fp(i+1) = fp(i+1) + L(i+1,k+1) * fcol(k+1);
    end
end
if isrow(f)
    fp = fp.';
end
end
