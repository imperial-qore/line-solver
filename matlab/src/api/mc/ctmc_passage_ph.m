function [alpha, S, s0, keep, atom] = ctmc_passage_ph(Q, pi0, target)
% [ALPHA, S, S0, KEEP, ATOM] = CTMC_PASSAGE_PH(Q, PI0, TARGET)
%
% Phase-type representation of the first passage time from the initial law PI0
% into the target state set TARGET, in the CTMC with generator Q.
%
% This is the primitive behind the whole ctmc_passage_* family. With A the
% complement of TARGET,
%
%     S     = Q(A,A)          sub-generator: the passage has not completed
%     s0    = -S*1  (= Q(A,TARGET)*1)   exit vector
%     alpha = PI0(A)          UNNORMALIZED, see below
%     atom  = sum(PI0(TARGET))
%
% so that L(s) = alpha (sI - S)^{-1} s0 + atom and F(t) = 1 - alpha exp(St) 1.
% Compare Eqs. 1-2 of P. G. Harrison and W. J. Knottenbelt, "Passage Time
% Distributions in Large Markov Chains", 2002, which write the same system as
% n linear equations with L_i = 1 on the target.
%
% ALPHA IS DELIBERATELY NOT NORMALIZED. Its mass is 1 - ATOM; the missing mass
% is the ATOM AT ZERO carried by initial states already in the target set. A
% caller that normalizes alpha and forgets the atom reports F(0) = 0 for a
% passage that has already completed with probability ATOM.
%
% PI0 may be empty, in which case the conditional stationary law on A is used.
% TARGET is 1-based here and 0-based in the Java, Python and C++ twins.
%
% KEEP maps the rows of S back to state indices of Q.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(Q,1);
if size(Q,2) ~= n
    line_error(mfilename, 'The generator must be square.');
end
target = unique(reshape(target, 1, []));
if isempty(target)
    line_error(mfilename, 'The target state set is empty: a first passage time into no state is undefined.');
end
if any(target < 1) || any(target > n)
    line_error(mfilename, 'A target state index is outside the state space.');
end
if max(abs(sum(Q,2))) > 1e-8 * max(1, max(abs(Q(:))))
    line_error(mfilename, 'Q is not an infinitesimal generator: its rows do not sum to zero. Pass it through ctmc_makeinfgen first.');
end

isTarget = false(1,n);
isTarget(target) = true;
keep = find(~isTarget);

S = Q(keep, keep);
s0 = -S * ones(length(keep),1);

if nargin < 2 || isempty(pi0)
    p = ctmc_solve(Q);
    p = reshape(p, 1, []);
    mass = sum(p(keep));
    if mass <= 0
        line_error(mfilename, 'The stationary law puts no mass outside the target set, so there is no passage to time.');
    end
    alpha = p(keep) / mass;
    atom = 0;
else
    pi0 = reshape(pi0, 1, []);
    if numel(pi0) ~= n
        line_error(mfilename, 'PI0 must be a distribution over the state space, one entry per state.');
    end
    alpha = pi0(keep);
    atom = sum(pi0(target));
end
end
