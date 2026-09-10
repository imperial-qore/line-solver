function b = moment_binomial_from_tail(t)
% b = moment_binomial_from_tail(t)
%
% Converts the survival (tail) probabilities t_m = P(N >= m) of a nonnegative
% integer random variable into its binomial moments,
%
%   b_j = E[nchoosek(N,j)] = sum_{m>=j} nchoosek(m-1,j-1) * t_m,  j >= 1
%
% with b_0 = t_0 = 1. Unlike every other edge of the house of moments this
% transform is UPPER triangular, so it consumes the whole tail: the result is
% exact only if the sequence covers the support, i.e. t_m = 0 beyond the last
% element supplied. This is the natural entry point for a closed queueing
% network, whose queue lengths are bounded by the population and whose joint
% survival probabilities are ratios of normalizing constants.
%
% Truncating the tail early yields a strict LOWER bound on every b_j, since all
% the coefficients and all the tail values are nonnegative. The bound is not
% inherited by the central moments downstream, whose conversion alternates in
% sign.
%
% Input:
%   t: vector of length n+1 holding t_0,...,t_n, i.e. t(i) is P(N >= i-1) and
%      t(1) = 1
%
% Output:
%   b: vector of length n+1 holding b_0,...,b_n, with the same orientation as t
%
% Example:
%   b = moment_binomial_from_tail([1, 1, 1, 1, 0]);  % N = 3 with probability 1
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

tcol = t(:);
n = length(tcol)-1;
b = moment_housematrix('binomial_from_tail', n) * tcol;
if isrow(t)
    b = b.';
end
end
