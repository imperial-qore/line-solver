function passage_firstpassage()
% PASSAGE_FIRSTPASSAGE
%
% First passage times in a Markov chain, and the exact cycle time along an
% overtake-free path of a closed tree-like product-form network.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% 1. The time for an M/M/1/K queue to fill from empty
% This is a first passage into a STATE SET, which no response-time getter can
% express: it is the chain reaching a marking, not a job finishing service.
K = 6; lambda = 1; mu = 1.5;
n = K+1;
Q = zeros(n);
for i = 1:n
    if i < n, Q(i,i+1) = lambda; end
    if i > 1, Q(i,i-1) = mu; end
end
Q = ctmc_makeinfgen(Q);

pi0 = zeros(1,n); pi0(1) = 1;      % start empty
target = n;                        % the full buffer

[mall, m] = ctmc_passage_moments(Q, pi0, target, 3);
fprintf('Time to fill an M/M/1/%d from empty (lambda=%g, mu=%g)\n', K, lambda, mu);
fprintf('  mean            = %.6f\n', m(1));
fprintf('  variance        = %.6f\n', m(2) - m(1)^2);
fprintf('  coeff. of var.  = %.6f\n', sqrt(m(2) - m(1)^2)/m(1));
fprintf('  from state K-1  = %.6f (one arrival away)\n', mall(n-1,1));

tset = linspace(0, 4*m(1), 400);
[F, f] = ctmc_passage_time(Q, pi0, target, tset);
fprintf('  P(fill <= mean) = %.6f\n', interp1(tset, F, m(1)));

% The mean hitting time from every state at once, the CTMC twin of
% dtmc_hitting_time.
h = ctmc_hitting_time(Q, target);
fprintf('  hitting times   = %s\n', mat2str(round(h(:)',3)));

%% 2. The same passage with a non-exponential sojourn (semi-Markov)
% The embedded chain is unchanged; only the holding-time law moves. Nothing in
% a generator can express this, which is why the semi-Markov route exists.
rate = -diag(Q)';
P = zeros(n);
for i = 1:n
    if rate(i) > 0
        P(i,:) = Q(i,:)/rate(i); P(i,i) = 0;
    else
        P(i,i) = 1;
    end
end
% Deterministic sojourns of the same mean: same embedded chain, tighter law.
hmom = zeros(n,3);
for i = 1:n
    d = 1/rate(i);
    hmom(i,:) = [d, d^2, d^3];
end
[~, mD] = smp_passage_moments(P, hmom, pi0, target, 3);
fprintf('\nSame embedded chain, DETERMINISTIC sojourns of equal mean:\n');
fprintf('  mean            = %.6f (unchanged, as it must be)\n', mD(1));
fprintf('  coeff. of var.  = %.6f (was %.6f)\n', ...
    sqrt(mD(2)-mD(1)^2)/mD(1), sqrt(m(2)-m(1)^2)/m(1));

%% 3. The cycle time of the tree-like network of Fig. 6 of the paper
mu6 = [3 5 4 6 2 1];
p12 = 0.2; p13 = 0.5; p14 = 0.3;
v = [1 p12 p13 p14 p12 p14];
N = 18;
paths = {[1 3], [1 2 5], [1 4 6]};
opt = struct('pathprob', [p13 p12 p14], 'nmom', 3);
tt = 0:0.25:40;
[fc, Fc, mom, out] = pfqn_cyclet_ofree(v, mu6, N, paths, tt, opt);
fprintf('\nTree network of Fig. 6, N = %d customers\n', N);
fprintf('  moments  = %.5f  %.4f  %.3f\n', mom(1), mom(2), mom(3));
fprintf('  paper    = 6.12717  53.3067  612.887\n');
fprintf('  routes   = %s\n', strjoin({out.method}, ', '));
fprintf('  P(cycle <= mean) = %.6f\n', interp1(tt, Fc, mom(1)));
end
