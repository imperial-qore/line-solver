function L = ctmc_passage_lst(Q, pi0, target, s)
% L = CTMC_PASSAGE_LST(Q, PI0, TARGET, S)
%
% Laplace-Stieltjes transform of the first passage time from PI0 into the
% target state set, evaluated at the (possibly complex) points S.
%
%     L(s) = alpha (sI - S)^{-1} s0 + atom
%
% which is Eqs. 1-2 of P. G. Harrison and W. J. Knottenbelt, "Passage Time
% Distributions in Large Markov Chains", 2002: one linear system per value of
% s, of the size of the non-target block.
%
% ONE SPARSE SOLVE PER s, NOT PER (s,t) PAIR. The saving over a dense matrix
% exponential is that the solves are sparse, so this route reaches chains a
% dense expm cannot hold. It is NOT a saving in the number of time points:
% every Abate-Whitt inverter places its nodes at s = beta/t, so a grid of T
% points costs T*|beta| solves. On a small chain ctmc_passage_time's default
% 'expm' route is faster.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[alpha, S, s0, ~, atom] = ctmc_passage_ph(Q, pi0, target);
nA = size(S,1);
L = zeros(size(s));
I = speye(nA);
for i = 1:numel(s)
    L(i) = alpha * ((s(i)*I - S) \ s0) + atom;
end
end
