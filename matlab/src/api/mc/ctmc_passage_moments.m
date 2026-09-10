function [mall, m] = ctmc_passage_moments(Q, pi0, target, nmax)
% [MALL, M] = CTMC_PASSAGE_MOMENTS(Q, PI0, TARGET, NMAX)
%
% Moments of order 1..NMAX of the first passage time into the target state
% set. MALL is (nstates x NMAX): row i is the moment vector for a passage
% STARTED IN STATE i, zero on target states and Inf where the target cannot be
% reached. M is the (1 x NMAX) moment vector for the initial law PI0.
%
% This is Eq. 3 of P. G. Harrison and W. J. Knottenbelt, "Passage Time
% Distributions in Large Markov Chains", 2002,
%
%     -q_ii M_i(n) = sum_{k not in target} q_ik M_k(n) + n M_i(n-1)
%
% i.e. (-S) M(n) = n M(n-1) with M(0) = 1, solved once per order: NMAX sparse
% solves and no transform inversion at all. The equivalent closed form is
% n! alpha (-S)^{-n} 1, which is NOT how it is evaluated here -- forming the
% inverse of the sub-generator destroys the sparsity the recursion preserves.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(nmax)
    nmax = 1;
end
[alpha, S, ~, keep, atom] = ctmc_passage_ph(Q, pi0, target);
n = size(Q,1);
nA = length(keep);

mall = zeros(n, nmax);
if nA == 0
    m = zeros(1, nmax);
    return
end

% A state that cannot reach the target has an infinite passage time; the
% sub-generator is singular on that block and a least-squares solve would
% return a finite number instead of saying so.
reach = local_reaches_target(Q, keep, n, target);

A = -S;
x = ones(nA,1);
X = zeros(nA, nmax);
warnState = warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
for k = 1:nmax
    x = A \ (k * x);
    X(:,k) = x;
end
warning(warnState);
X(~reach, :) = Inf;

mall(keep, :) = X;
m = zeros(1, nmax);
finite = reach(:)' & (alpha(:)' ~= 0);
if any(~reach(:)' & alpha(:)' > 0)
    m(:) = Inf;
else
    for k = 1:nmax
        m(k) = alpha(finite) * X(finite, k);
    end
end
% States already in the target contribute a zero passage time, so the atom
% lowers no moment: it is recorded for the caller rather than folded in.
if atom >= 1
    m(:) = 0;
end
end

function reach = local_reaches_target(Q, keep, n, target)
% Backward reachability closure over the transition graph.
adj = (Q ~= 0);
adj(1:n+1:end) = false;
seen = false(1,n);
seen(target) = true;
frontier = target;
while ~isempty(frontier)
    pred = find(any(adj(:, frontier), 2))';
    pred = pred(~seen(pred));
    seen(pred) = true;
    frontier = pred;
end
reach = seen(keep)';
end
