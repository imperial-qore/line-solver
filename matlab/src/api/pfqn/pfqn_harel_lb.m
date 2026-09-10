%{
%{
 % @file pfqn_harel_lb.m
 % @brief Harel-Namn-Sturm throughput lower bound of a closed network.
%}
%}

%{
%{
 % @brief Harel-Namn-Sturm throughput lower bound of a closed network.
 % @fn pfqn_harel_lb(rho, N, Z)
 % @param rho Relative utilizations (k x 1), all strictly positive.
 % @param N Closed population, at least 1.
 % @param Z Think time; must be zero.
 % @return LB Throughput lower bound at population N.
%}
%}
function LB = pfqn_harel_lb(rho, N, Z)
% LB = PFQN_HAREL_LB(RHO, N, Z)
%
% Lower bound alone of PFQN_HAREL_BOUNDS,
%
%   LB = N / (A_1 + (N-1) (A_N/A_1)^{1/(N-1)}),   A_i = sum_j rho_j^i,
%
% from Harel, Namn and Sturm, "Simple bounds for closed queueing networks"
% (Queueing Systems 31, 1999). A nonzero think time is refused.
%
% See also PFQN_HAREL_BOUNDS, PFQN_HAREL_UB.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(Z)
    Z = 0;
end
if Z ~= 0
    line_error(mfilename, ...
        'pfqn_harel_lb is only valid for networks with zero think time; the provided think time is nonzero.');
end
if N < 1
    line_error(mfilename, 'pfqn_harel_lb: the population must be at least 1.');
end
rho = rho(:);
if isempty(rho)
    line_error(mfilename, 'pfqn_harel_lb: the loading vector must have at least one element.');
end
if any(rho <= 0)
    line_error(mfilename, 'pfqn_harel_lb: all loading factors must be positive.');
end

A1 = sum(rho);
if N == 1
    LB = 1 / A1;
    return
end
AN = sum(rho .^ N);
LB = N / (A1 + (N - 1) * (AN / A1) ^ (1 / (N - 1)));
end
