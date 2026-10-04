%{
%{
 % @file pfqn_harel_ub.m
 % @brief Harel-Namn-Sturm throughput upper bound of a closed network.
%}
%}

%{
%{
 % @brief Harel-Namn-Sturm throughput upper bound of a closed network.
 % @fn pfqn_harel_ub(rho, N, n, Z)
 % @param rho Relative utilizations (k x 1), all strictly positive.
 % @param N Closed population, at least 1.
 % @param n Extrapolation point, 2 <= n <= min(N,7).
 % @param Z Think time; must be zero.
 % @return UB Throughput upper bound at population N.
%}
%}
function UB = pfqn_harel_ub(rho, N, n, Z)
% UB = PFQN_HAREL_UB(RHO, N, N0, Z)
%
% Upper bound alone of PFQN_HAREL_BOUNDS, extrapolated from the EXACT
% throughput TH(n) = G(n-1)/G(n) at the small population n,
%
%   UB(n) = N / (A_1 + ((N-1)/(n-1)) (n/TH(n) - A_1)),   A_i = sum_j rho_j^i,
%
% from Harel, Namn and Sturm, "Simple bounds for closed queueing networks"
% (Queueing Systems 31, 1999). G is evaluated by the Newton-Girard recurrence
% n G(n) = sum_{i=1..n} A_i G(n-i); the n <= 7 ceiling is kept from the
% reference implementation. A nonzero think time is refused.
%
% See also PFQN_HAREL_BOUNDS, PFQN_HAREL_LB.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(Z)
    Z = 0;
end
if Z ~= 0
    line_error(mfilename, ...
        'pfqn_harel_ub is only valid for networks with zero think time; the provided think time is nonzero.');
end
if N < 1
    line_error(mfilename, 'pfqn_harel_ub: the population must be at least 1.');
end
if n < 2
    line_error(mfilename, 'pfqn_harel_ub: the extrapolation point must be at least 2.');
end
if n > N
    line_error(mfilename, 'pfqn_harel_ub: the extrapolation point cannot exceed N.');
end
if n > 7
    line_error(mfilename, 'pfqn_harel_ub: the extrapolation point cannot exceed 7.');
end
rho = rho(:);
if isempty(rho)
    line_error(mfilename, 'pfqn_harel_ub: the loading vector must have at least one element.');
end
if any(rho <= 0)
    line_error(mfilename, 'pfqn_harel_ub: all loading factors must be positive.');
end

A = zeros(1, n);
for i = 1:n
    A(i) = sum(rho .^ i);
end
G = zeros(1, n + 1);
G(1) = 1;
for m = 1:n
    acc = 0;
    for i = 1:m
        acc = acc + A(i) * G(m - i + 1);
    end
    G(m + 1) = acc / m;
end
if G(n + 1) == 0
    line_error(mfilename, 'pfqn_harel_ub: the normalizing constant vanishes.');
end
THn = G(n) / G(n + 1);
den = A(1) + ((N - 1) / (n - 1)) * (n / THn - A(1));
if den == 0
    line_error(mfilename, 'pfqn_harel_ub: the upper-bound denominator vanishes.');
end
UB = N / den;
end
