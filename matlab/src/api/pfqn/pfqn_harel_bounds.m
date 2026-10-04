%{
%{
 % @file pfqn_harel_bounds.m
 % @brief Harel-Namn-Sturm throughput bounds for a single-class closed network.
%}
%}

%{
%{
 % @brief Harel-Namn-Sturm throughput bounds for a single-class closed network.
 % @fn pfqn_harel_bounds(rho, N, Z, maxUB)
 % @param rho Relative utilizations (k x 1), all strictly positive.
 % @param N Closed population, at least 1.
 % @param Z Think time; must be zero.
 % @param maxUB Largest extrapolation point; defaults to min(N,7).
 % @return LB Throughput lower bound at population N.
 % @return UB Upper bounds UB(n) for n = 2..maxUB (UB(1) unset).
 % @return TH Exact throughput TH(n) at population n, n = 1..maxUB.
%}
%}
function [LB, UB, TH] = pfqn_harel_bounds(rho, N, Z, maxUB)
% [LB, UB, TH] = PFQN_HAREL_BOUNDS(RHO, N, Z, MAXUB)
%
% Harel-Namn-Sturm throughput bounds of a single-class closed network of
% load-independent queues with relative utilizations rho.
%
% These are the SHARP bounds of Harel, Namn and Sturm, "Simple bounds for
% closed queueing networks" (Queueing Systems 31, 1999), distinct from the
% 'sb' family already in SOLVER_BA: 'sb' uses only the first three power sums
% in closed form, whereas this family evaluates the normalizing constant
% exactly at small populations and extrapolates from it. Both cite the same
% paper; they are different results in it and neither subsumes the other.
%
% With the power sums A_i = sum_j rho_j^i,
%
%   G(n)  = h_n(rho), the complete homogeneous symmetric polynomial,
%   TH(n) = G(n-1)/G(n),                 the exact throughput at population n,
%   LB    = N / (A_1 + (N-1) (A_N/A_1)^{1/(N-1)}),
%   UB(n) = N / (A_1 + ((N-1)/(n-1)) (n/TH(n) - A_1)),        2 <= n <= N.
%
% G(n) IS the normalizing constant of the closed load-independent network at
% population n, so it must equal PFQN_CA on the same demands and TH(n) must
% equal the exact PFQN_MVA throughput at population n; that identity is the
% oracle the implementation is tested against. G is evaluated by the
% Newton-Girard recurrence n G(n) = sum_{i=1..n} A_i G(n-i) rather than from
% expanded polynomials in the power sums; the n <= 7 ceiling on the
% extrapolation point is kept from the reference implementation.
%
% A nonzero think time is REFUSED rather than folded in: the bounds are
% derived for a network with no terminal population, so silently dropping Z
% would return a bound that does not bound.
%
% See also PFQN_HAREL_LB, PFQN_HAREL_UB, PFQN_CA, PFQN_MVA, SOLVER_BA_ANALYZER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(Z)
    Z = 0;
end
if nargin < 4
    maxUB = [];
end
harel_reject_thinktime(Z, 'pfqn_harel_bounds');
if N < 1
    line_error(mfilename, 'pfqn_harel_bounds: the population must be at least 1.');
end
rho = rho(:);
harel_check_rho(rho);

if isempty(maxUB) || maxUB <= 0
    maxUB = min(N, 7);
end
if maxUB > 7
    line_error(mfilename, 'pfqn_harel_bounds: upper bounds are available only for n <= 7.');
end
if maxUB > N
    line_error(mfilename, 'pfqn_harel_bounds: the extrapolation point cannot exceed N.');
end

% The lower bound reads A up to N, the upper bounds only up to maxUB.
% A(i) is the power sum of order i; G(n+1) is G(n), with G(1) = G(0) = 1.
A = harel_power_sums(rho, max(N, maxUB));
LB = harel_lower_bound(A, N);

G = harel_G(A, maxUB);
TH = zeros(1, maxUB);
UB = zeros(1, maxUB);
for n = 1:maxUB
    if G(n + 1) == 0
        line_error(mfilename, 'pfqn_harel_bounds: the normalizing constant vanishes.');
    end
    TH(n) = G(n) / G(n + 1);
end
for n = 2:maxUB
    UB(n) = harel_upper_from_th(A(1), N, n, TH(n));
end
end

% ---- power sums A_i = sum_j rho_j^i, i = 1..maxPower ----------------------
function A = harel_power_sums(rho, maxPower)
A = zeros(1, maxPower);
for i = 1:maxPower
    A(i) = sum(rho .^ i);
end
end

% ---- G(0..n) by the Newton-Girard recurrence n G(n) = sum_i A_i G(n-i) ----
function G = harel_G(A, n)
if numel(A) < n
    line_error('pfqn_harel_bounds', 'too few power sums for the requested population.');
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
end

% ---- LB = N / (A1 + (N-1) (A_N/A_1)^{1/(N-1)}) ----------------------------
function LB = harel_lower_bound(A, N)
if N == 1
    LB = 1 / A(1);
    return
end
LB = N / (A(1) + (N - 1) * (A(N) / A(1)) ^ (1 / (N - 1)));
end

% ---- UB(n) = N / (A1 + ((N-1)/(n-1)) (n/TH(n) - A1)) ----------------------
function UB = harel_upper_from_th(A1, N, n, THn)
if THn == 0
    line_error('pfqn_harel_bounds', 'the throughput at the extrapolation point is zero.');
end
den = A1 + ((N - 1) / (n - 1)) * (n / THn - A1);
if den == 0
    line_error('pfqn_harel_bounds', 'the upper-bound denominator vanishes.');
end
UB = N / den;
end

% ---- the reference refuses a nonzero think time rather than folding it in --
function harel_reject_thinktime(Z, who)
if Z ~= 0
    line_error(who, ...
        '%s is only valid for networks with zero think time; the provided think time is nonzero.', who);
end
end

% ---- shared input screening of the loading vector -------------------------
function harel_check_rho(rho)
if isempty(rho)
    line_error('pfqn_harel_bounds', 'the loading vector must have at least one element.');
end
if any(rho <= 0)
    line_error('pfqn_harel_bounds', 'all loading factors must be positive.');
end
end
