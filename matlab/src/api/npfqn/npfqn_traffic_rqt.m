function [lambda,Gamma,alpha] = npfqn_traffic_rqt(lambda0,Gamma0,alpha0,F)
% [LAMBDA,GAMMA,ALPHA] = NPFQN_TRAFFIC_RQT(LAMBDA0,GAMMA0,ALPHA0,F)
%
% Effective arrival process perceived at each node of a single-class open
% queueing network under the Robust Queueing Theory (RQT) calculus. Solves the
% network characterization of Theorem 10 (Theorem 7 when all tail coefficients
% agree), which composes three operators: passage through a queue with
% adversarial servers leaves the uncertainty set unchanged (robust Burke,
% Theorem 4), superposition merges sets by Theorem 5, and thinning by a fraction
% f scales the rate by f and the variability by f^(-1/alpha) (Theorem 6). The
% resulting equations are
%
%   lambda_j = lambda0_j + sum_i lambda_i f_ij,
%   Gamma_j  = (1/lambda_j) [ 1{a0_j=ab_j} (lambda0_j Gamma0_j)^(p_j)
%                             + sum_i 1{ab_i=ab_j} (lambda_i Gamma_i)^(p_i) f_ij ]^(1/p_j),
%
% with p_j = ab_j/(ab_j-1) and ab_j = min over the streams feeding j of their
% tail coefficients: the heaviest tail upstream dominates. Both are solved
% exactly rather than iteratively. The rate equations are the usual traffic
% equations, and in the variables z_j = (lambda_j Gamma_j)^(p_j) the variability
% equations are linear as well, so each is one linear system; ab is obtained by
% propagating the minimum to a fixed point.
%
% Inputs:
%   LAMBDA0 - (J x 1) external arrival rate at each node, 0 where there is none
%   GAMMA0  - (J x 1) variability parameter of each external arrival process,
%             which for a renewal stream is the interarrival standard deviation
%   ALPHA0  - (J x 1) tail coefficient in (1,2] of each external arrival process
%   F       - (J x J) routing probability matrix, F(i,j) = fraction of the jobs
%             leaving node i that go to node j (row sums <= 1)
%
% Returns:
%   LAMBDA  - (J x 1) effective arrival rate at each node
%   GAMMA   - (J x 1) effective variability parameter at each node
%   ALPHA   - (J x 1) effective tail coefficient at each node
%
% Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
% Operations Research 63(3), 676-700, Theorems 4-7 and 10.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lambda0 = lambda0(:); Gamma0 = Gamma0(:); alpha0 = alpha0(:);
J = numel(lambda0);
F = full(F);

% traffic equations
lambda = (eye(J) - F') \ lambda0;
lambda(abs(lambda) < eps) = 0;

% effective tail coefficient: the minimum propagated along the routing graph
alpha = repmat(2, J, 1);
hasExt = lambda0 > 0;
alpha(hasExt) = alpha0(hasExt);
alpha(~hasExt) = Inf;
for it = 1:J
    prev = alpha;
    for j = 1:J
        for i = 1:J
            if F(i,j) > 0 && lambda(i) > 0
                alpha(j) = min(alpha(j), alpha(i));
            end
        end
    end
    if isequal(prev, alpha)
        break
    end
end
alpha(~isfinite(alpha)) = 2; % an unreachable node keeps the light-tailed default

% variability equations, linear in z_j = (lambda_j Gamma_j)^(p_j)
p = alpha ./ (alpha - 1);
z0 = zeros(J,1);
for j = 1:J
    if lambda0(j) > 0 && abs(alpha0(j) - alpha(j)) < 1e-12
        z0(j) = (lambda0(j) * Gamma0(j))^p(j);
    end
end
A = zeros(J,J);
for i = 1:J
    for j = 1:J
        if F(i,j) > 0 && abs(alpha(i) - alpha(j)) < 1e-12
            A(i,j) = F(i,j);
        end
    end
end
z = (eye(J) - A') \ z0;
z(z < 0) = 0;

Gamma = zeros(J,1);
for j = 1:J
    if lambda(j) > 0
        Gamma(j) = z(j)^(1/p(j)) / lambda(j);
    end
end
end
