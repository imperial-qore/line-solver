function test_qrf_noblo_bethe
% Regression test for QRF_NOBLO_BETHE, the tree-reweighted (Bethe) arm.
%
% No LINE dependencies and no LP solver: the oracles are the exact
% product-form law of a two-station cycle, the polytope invariants read off
% the returned tensor, and the native-Python reference on a three-station
% instance where the polytope has slack.
%
% WHAT IS ASSERTED, AND WHY IT IS NOT A TABLE OF DIGITS
% -----------------------------------------------------
% QRF_NOBLO_BETHE minimises lambda*sum_{i~=j} I(n_i;n_j) - sum_i H(n_i) with
% lambda = 1/M, the negative of a tree-reweighted entropy at the uniform point
% rho_ij = 2/M of the spanning tree polytope of K_M. That weight is the largest
% uniform one at which the program is CONVEX, and what convexity buys is
% well-posedness, not accuracy: every local optimum is global, so the answer is
% a property of the model rather than of the start point.
%
% Digits are kept out of the loose-polytope row on purpose. Restoring the n = 0
% cells -- the range the AMPL source states and the coded mmi() does not use --
% brings the structurally zero entries inside the sum, where they contribute
% 0*log(LOGTOL) = 0 to the objective VALUE but log(LOGTOL) ~ -13.8 to the
% GRADIENT. The value is insensitive to LOGTOL, the descent direction is not,
% and fmincon here differs from the ports' SLSQP and Frank-Wolfe anyway, so the
% cross-codebase row is asserted at 1e-3 and the pinned rows at 1e-6.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% 1. THE PINNED POLYTOPE. At M == 2 the pairwise joint of a closed chain is
% fixed by the marginal, so the equality rows leave no freedom and EVERY
% objective over this polytope must report the exact product-form answer. This
% is the row that catches a broken objective as a wrong number rather than as a
% slow solve.
pinned = {struct('mu',[1.0 2.0],'N',2), ...
          struct('mu',[1.0 1.5],'N',3), ...
          struct('mu',[4.0 1.0],'N',3)};
for t = 1:numel(pinned)
    mu2s = pinned{t}.mu;
    N = pinned{t}.N;
    M = 2;
    rt = [0 1; 1 0];
    K = [1 1];
    mu = zeros(M,1,1); v = zeros(M,1,1);
    mu(1,1,1) = mu2s(1);
    mu(2,1,1) = mu2s(2);

    [UN,QN,p2opt] = qrf_noblo_bethe(M, 1, K, N, mu, v, rt);
    [UNex,QNex] = local_exact_cycle(mu2s(1), mu2s(2), N);

    assert(max(abs(UN(:)' - UNex)) < 1e-6, ...
        'pinned case %d: U = %s, exact %s', t, mat2str(UN,8), mat2str(UNex,8));
    assert(max(abs(QN(:)' - QNex)) < 1e-6, ...
        'pinned case %d: Q = %s, exact %s', t, mat2str(QN,8), mat2str(QNex,8));
    assert(abs(sum(QN) - N) < 1e-6, ...
        'pinned case %d: population %g, expected %d', t, sum(QN), N);
    local_check_feasible(p2opt, M, N, K, 1, sprintf('pinned case %d', t));
end

%% 2. A POLYTOPE WITH SLACK. Three exponential stations at N = 3: the free
% dimension is 3, so the objective -- not the equality rows -- decides the
% answer, which is what makes this the row that would move if the population
% loops or the sign of the entropy term were wrong.
M = 3; N = 3;
K = [1 1 1];
rates = [1.0 1.5 2.0];
rt = zeros(M);
for i = 1:M
    rt(i, mod(i,M)+1) = 1;
end
mu = zeros(M,1,1); v = zeros(M,1,1);
for i = 1:M
    mu(i,1,1) = rates(i);
end

[UN,QN,p2opt] = qrf_noblo_bethe(M, 1, K, N, mu, v, rt);

% Native Python `qrf_noblo_bethe(3,1,[1 1 1],3,mu,v,rt)` on THIS instance, and
% the JAR agrees to the digits printed. MATLAB solves with fmincon and the
% ports with SLSQP / conditional gradient, so the tolerance is 1e-3 rather than
% solver precision. The glpsol range of U1 over the same polytope is
% [0.77593985, 0.82438580], so the value below is a strict interior point of it
% -- the method is an estimator here, not a bound.
UNref = [0.79768565 0.53179043 0.39884283];
assert(max(abs(UN(:)' - UNref)) < 1e-3, ...
    'slack case: U = %s, native-Python reference %s', mat2str(UN,8), mat2str(UNref,8));
assert(abs(sum(QN) - N) < 1e-6, 'slack case: population %g, expected %d', sum(QN), N);
assert(all(UN > 1e-6) && all(UN < 1 - 1e-6), ...
    'slack case: U = %s is on the boundary, which reads as a stalled solve', mat2str(UN,8));
% ... and the point it reports is a point of the polytope, not merely a vector.
local_check_feasible(p2opt, M, N, K, 1, 'slack case');

fprintf('test_qrf_noblo_bethe: all checks passed\n');
end

% ---------------------------------------------------------------------------

function [UN,QN] = local_exact_cycle(mu1, mu2, N)
% The exact law of the two-station cyclic M/M/1//N network.
x = [1/mu1, 1/mu2];
w = zeros(1,N+1);
for n1 = 0:N
    w(n1+1) = x(1)^n1 * x(2)^(N-n1);
end
p = w / sum(w);
UN = [sum(p(2:end)), sum(p(1:end-1))];
QN = [sum((0:N).*p), sum((N:-1:0).*p)];
end

function local_check_feasible(p2, M, N, K, MR, label)
% The polytope invariants that can be read off the returned tensor alone.
%
% ONE, SYMMETRY, ZERO2, ZERO3 and MARGINALS are the blocks that define the
% object as a locally consistent pairwise law; THM2 and COR1 are the two
% population identities. A solve that returned a vector rather than a point of
% the polytope fails here whatever its objective value.
tol = 1e-6;

% ONE: sum over the diagonal cells of each station is 1
for j = 1:M
    s = 0;
    for nj = 1+(0:N), for k = 1:K(j), for m = 1:MR
        s = s + p2(j,nj,k,j,nj,k,m);
    end, end, end
    assert(abs(s-1) < tol, '%s: ONE violated at station %d by %g', label, j, abs(s-1));
end

% SYMMETRY, ZERO2, ZERO3
for j = 1:M, for nj = 1+(0:N), for k = 1:K(j)
    for i = 1:M, for ni = 1+(0:N), for h = 1:K(i), for m = 1:MR
        assert(abs(p2(i,ni,h,j,nj,k,m) - p2(j,nj,k,i,ni,h,m)) < tol, ...
            '%s: SYMMETRY violated', label);
        if i==j && ni~=nj
            assert(abs(p2(j,nj,k,i,ni,h,m)) < tol, '%s: ZERO2 violated', label);
        end
        if i~=j && (nj-1)+(ni-1) > N
            assert(abs(p2(j,nj,k,i,ni,h,m)) < tol, '%s: ZERO3 violated', label);
        end
    end, end, end, end
end, end, end

% MARGINALS: the diagonal cell is the marginal of every pairing
for j = 1:M, for k = 1:K(j), for nj = 1+(0:N), for i = 1:M, for m = 1:MR
    if i ~= j
        s = 0;
        for ni = 1+(0:N), for h = 1:K(i)
            s = s + p2(j,nj,k,i,ni,h,m);
        end, end
        d = abs(p2(j,nj,k,j,nj,k,m) - s);
        assert(d < tol, '%s: MARGINALS violated by %g', label, d);
    end
end, end, end, end, end

% THM2: the conditional population is N
for j = 1:M, for k = 1:K(j), for nj = 1+(0:N), for m = 1:MR
    s = 0;
    for i = 1:M, for ni = 1+(1:N), for ki = 1:K(i)
        s = s + (ni-1)*p2(j,nj,k,i,ni,ki,m);
    end, end, end
    d = abs(s - N*p2(j,nj,k,j,nj,k,m));
    assert(d < tol, '%s: THM2 violated by %g', label, d);
end, end, end, end

% COR1: the second moment of the total population is N^2
s = 0;
for m = 1:MR, for i = 1:M, for j = 1:M
    for nj = 1+(1:N), for ni = 1+(1:N), for ki = 1:K(i), for kj = 1:K(j)
        s = s + (ni-1)*(nj-1)*p2(j,nj,kj,i,ni,ki,m);
    end, end, end, end
end, end, end
assert(abs(s - N^2) < 1e-5, '%s: COR1 violated by %g', label, abs(s-N^2));

% Every entry is a probability
assert(min(p2(:)) > -tol && max(p2(:)) < 1+tol, '%s: an entry left [0,1]', label);
end
