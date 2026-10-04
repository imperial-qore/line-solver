function en = spn_rec_enabled(mdds, g, mde, nplacelevels)
% EN = SPN_REC_ENABLED(MDDS, G_L, MDE, NPLACELEVELS)
% Enabling-degree distribution of one mode of a product-form stochastic Petri
% net, by the masked MDD-rec recursion.
%
% The enabling degree of a mode in marking m is
%
%   e(m) = min_{l : I_l > 0} floor(m_l / I_l),
%
% zero when any inhibitor threshold is met. P(e >= k) is therefore the mass of
% the marking subset in which EVERY input level holds at least k*I_l tokens and
% no inhibitor fires, which is a per-level restriction and so exactly what
% MDD_REC_MASKED computes: the paper's second modified recurrence is the same
% walk under a different mask, not a second algorithm.
%
% -- Input
% MDDS         : MDD.toStruct of the reachable set built by SPN_MDD
% G_L          : 1 x K cell of per-level product-form factors
% MDE          : one entry of INFO.modes as returned by SPN_MDD
% NPLACELEVELS : how many leading levels are place levels
% -- Output
% EN : struct with fields
%        ge        - ge(k+1) is the mass of {e >= k}; ge(1) is the whole set
%        eq        - eq(k+1) is the mass of {e == k}, i.e. ge(k+1)-ge(k+2)
%        maxDegree - the largest enabling degree the place bounds permit, E_j
%
% The masses are UNNORMALISED, as in the paper; divide by G from MDD_REC for
% probabilities. SPN_METRICS does that and turns them into the transition
% measures.
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 5.3.
%
% See also MDD_REC_MASKED, SPN_METRICS, SPN_MDD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = mdds.K;
L = nplacelevels;
if L > K
    line_error(mfilename, 'more place levels than diagram levels');
end

% E_j: the enabling degree cannot exceed what the tightest input place bound
% allows. A mode with no input place has no bound and is refused rather than
% silently truncated, matching SPN_MDD's own refusal.
emax = 0;
hasInput = false;
for l = 1:L
    if ~(mde.enab(l) > 0), continue; end
    cap = floor((mdds.domain(l) - 1) / mde.enab(l));
    if hasInput, emax = min(emax, cap); else, emax = cap; end
    hasInput = true;
end
if ~hasInput
    line_error(mfilename, ['the mode consumes from no place, so its enabling degree ' ...
        'is unbounded']);
end

ge = zeros(1, emax + 2);
for k = 0:emax
    mask = cell(1, K);
    for j = 1:K, mask{j} = true(1, mdds.domain(j)); end
    for l = 1:L
        v = (0:mdds.domain(l) - 1);
        % k = 0 asks only that the marking exist, so the inhibitor test belongs
        % to k >= 1: e = 0 covers the inhibited markings too.
        drop = v < mde.enab(l) * k;
        if k > 0, drop = drop | (v >= mde.inhib(l)); end
        mask{l}(drop) = false;
    end
    ge(k + 1) = mdd_rec_masked(mdds, g, mask);
end
eq = zeros(1, emax + 2);
eq(1:emax + 1) = ge(1:emax + 1) - ge(2:emax + 2);

en = struct('ge', ge, 'eq', eq, 'maxDegree', emax);
end
