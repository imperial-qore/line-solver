function t = moment_tail_from_binomial(b)
% t = moment_tail_from_binomial(b)
%
% Converts the binomial moments of a nonnegative integer random variable into
% its survival (tail) probabilities,
%
%   t_m = sum_{j>=m} (-1)^(j-m) * nchoosek(j-1,m-1) * b_j,  m >= 1
%
% with t_0 = 1. Inverse of moment_binomial_from_tail. The inversion is exact on
% the finite box supplied, the matrix being unit upper triangular, but it
% reconstructs the true tail only if the binomial moments were themselves those
% of a law supported on 0,...,n.
%
% Input:
%   b: vector of length n+1 holding b_0,...,b_n
%
% Output:
%   t: vector of length n+1 holding t_0,...,t_n, with the same orientation as b
%
% Example:
%   t = moment_tail_from_binomial(moment_binomial_from_tail([1, 1, 1, 0]));
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

bcol = b(:);
n = length(bcol)-1;
t = moment_housematrix('tail_from_binomial', n) * bcol;
if isrow(b)
    t = t.';
end
end
