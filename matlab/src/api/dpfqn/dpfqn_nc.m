%{
%{
 % @file pfqn_ncdt.m
 % @brief Normalizing constants of a discrete-time closed cycle of Bernoulli
 %        servers with state independent service probabilities.
%}
%}

function [lG, G, G1] = pfqn_ncdt(p, N)
%{
%{
 % @brief Buzen-style recursions for the discrete-time closed cycle.
 % @fn pfqn_ncdt(p, N)
 % @param p Row vector of per-slot service completion probabilities, p(j) in (0,1).
 % @param N Number of customers cycling in the J nodes.
 % @return lG Logarithm of the time-stationary normalizing constant G(N,J).
 % @return G Time-stationary normalizing constant G(N,J).
 % @return G1 Column vector of arrival normalizing constants, G1(k+1)=G_1(k,J).
%}
%}
% PFQN_NCDT computes the two families of normalizing constants of a closed
% cycle of J state independent Bernoulli servers observed on a discrete time
% scale. With q(j)=1-p(j) the time-stationary law of the queue length vector is
%
%   pi(n_1,...,n_J) = prod_j (q_j/p_j)^n_j (1/q_j)^{1{n_j>0}} / G(N,J),
%
% i.e. the discrete-time product form of Daduna (2001), corollary 3.4, whose
% extra factor (1/q_j)^{1{n_j>0}} on the busy nodes is what distinguishes it
% from the continuous-time Gordon-Newell form. G(N,J) obeys the three-term
% recursion of proposition 3.18
%
%   G(k,j) = G(k,j-1) + (q_j/p_j) G(k-1,j) + G(k-1,j-1),
%
% with G(0,j)=1 and G(k,0)=0 for k>=1. Unlike the continuous-time convolution
% algorithm this recursion is not invariant to the numbering of the nodes.
%
% G1 collects the arrival normalizing constants of proposition 3.17(b), the
% constants of the distribution seen by a customer joining a node, in which
% the arrival node carries the weight (q/p)^n with no 1/q factor. They obey
% proposition 3.19,
%
%   G1(k,J) = G(k-1,J-1) + (q_J/p_J) G1(k-1,J),   k >= 3,
%
% with G1(1,J)=1 and G1(2,J)=q_1/p_1 + sum_{j>=2} 1/p_j. By lemma 7.3 of the
% same reference the arrival constant is the same at every node, so one family
% suffices. Three consequences are what the discrete-time analyzer consumes:
%
%   throughput per slot   X   = G1(N,J) / G(N,J)   (equal at every node),
%   utilization           U_j = X / p_j,
%   tail probability      P(X_j >= k) = (q_j/p_j)^k (1/q_j) G1(N-k+1,J)/G(N,J).
%
% The last identity is corollary 3.20(a) of the reference with its index
% corrected: as printed there the right-hand side evaluates to P(X_j >= k+1).
% Corollary 3.20(c), which transfers the tail from node 1 to node j, holds for
% k >= 1 only; at k = 0 both tails are 1 while the stated ratio is q_1/q_j.
%
% G and G1 are returned in a common scale, which is unity unless the recursion
% had to be rescaled to keep (q/p)^N inside double range. Ratios such as
% G1(k)/G are exact in either case, and lG is always the true log constant.
%
% Examples:
%   [lG, G, G1] = pfqn_ncdt([0.5 0.25 0.7], 5);
%   G1(end)/G                    % throughput per slot, 0.2483
%
% Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
% Springer 2001, corollary 3.4, propositions 3.17-3.19 and corollary 3.20.
%
% See also PFQN_NCDTLD, QSYS_BERNOULLI1, PFQN_CA
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isnumeric(p) || ~isreal(p) || isempty(p) || any(~isfinite(p(:))) || any(p(:) <= 0) || any(p(:) >= 1)
    line_error(mfilename, 'service probabilities must be real and in the open interval (0,1)');
end
if ~isnumeric(N) || ~isscalar(N) || ~isreal(N) || N < 0 || floor(N) ~= N
    line_error(mfilename, 'N must be a non-negative integer');
end
p = p(:)';
J = numel(p);
x = (1 - p) ./ p;                 % per-node geometric ratio q_j/p_j

% Proposition 3.18. The recursion is linear and homogeneous in the whole
% table, so rescaling every entry at once preserves it; that is what keeps
% (q/p)^N from overflowing for a slow node and a large population. G1 is
% rescaled in step so that the two families stay in one common scale.
Gt = zeros(N+1, J+1);
Gt(1, :) = 1;                     % G(0,j) = 1, including j = 0
G1 = zeros(N+1, 1);               % G1(1) = G_1(0,J) = 0, no composition of -1
lscale = 0;
for k = 1:N
    for j = 1:J
        Gt(k+1, j+1) = Gt(k+1, j) + x(j) * Gt(k, j+1) + Gt(k, j);
    end
    switch k
        case 1
            G1(2) = Gt(1, 1);                                  % G_1(1,J) = 1
        case 2
            G1(3) = (x(1) + sum(1 ./ p(2:end))) * Gt(1, 1);
        otherwise
            G1(k+1) = Gt(k, J) + x(J) * G1(k);
    end
    mx = max(Gt(k+1, :));
    if mx > 1e250
        Gt = Gt / mx;
        G1 = G1 / mx;
        lscale = lscale + log(mx);
    end
end
G = Gt(N+1, J+1);
lG = log(G) + lscale;
end
