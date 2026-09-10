function f = moment_factorial_from_binomial(b)
% f = moment_factorial_from_binomial(b)
%
% Converts the binomial moments b_n = E[nchoosek(N,n)] of a discrete random
% variable N into the factorial moments f_n = E[N(N-1)...(N-n+1)] via the
% one-to-one correspondence
%
%   f_n = n! * b_n
%
% Input:
%   b: vector of length n+1 holding b_0,...,b_n, i.e. b(i) is the moment of
%      order i-1 and b(1) = b_0 = 1
%
% Output:
%   f: vector of length n+1 holding f_0,...,f_n, with the same orientation
%      as b
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003, Section 4.
%
% Example:
% f = moment_factorial_from_binomial([1,2,2,4/3])

bcol = b(:);
n = length(bcol)-1;
f = zeros(n+1,1);
for i = 0:n
    f(i+1) = factorial(i) * bcol(i+1);
end
if isrow(b)
    f = f.';
end
end
