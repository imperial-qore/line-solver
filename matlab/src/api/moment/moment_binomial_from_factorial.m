function b = moment_binomial_from_factorial(f)
% b = moment_binomial_from_factorial(f)
%
% Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete
% random variable N into the binomial moments b_n = E[nchoosek(N,n)] via the
% one-to-one correspondence
%
%   b_n = f_n / n!
%
% Input:
%   f: vector of length n+1 holding f_0,...,f_n, i.e. f(i) is the moment of
%      order i-1 and f(1) = f_0 = 1
%
% Output:
%   b: vector of length n+1 holding b_0,...,b_n, with the same orientation
%      as f
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% b = moment_binomial_from_factorial([1,2,4,8])

fcol = f(:);
n = length(fcol)-1;
b = zeros(n+1,1);
for i = 0:n
    b(i+1) = fcol(i+1) / factorial(i);
end
if isrow(f)
    b = b.';
end
end
