function mass = mdd_rec_marginal(mdds, g, l)
% MASS = MDD_REC_MARGINAL(MDDS, G_L, L)
% Unnormalised masses of {s in S : s_L = k}, one per local value k of level L.
%
% Divided by the normalising constant these are P(m_l = k) of Sec. 5.3 of the
% MDD-rec paper: the mean occupancy of a level is sum_k k * P(m_l = k), and its
% utilization 1 - P(m_l = 0).
%
% -- Input
% MDDS : MDD.toStruct of the reachable set
% G_L  : 1 x K cell, G_L{l}(v+1) is g_l(v)
% L    : level index, 1..K
% -- Output
% MASS : 1 x domain(L) unnormalised masses, MASS(k+1) for local value k
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 5.3.
%
% See also MDD_REC, MDD_REC_MASKED, SPN_METRICS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = mdds.K;
if l < 1 || l > K
    line_error(mfilename, 'level index is out of range');
end
d = mdds.domain(l);
mass = zeros(1, d);
for k = 1:d
    mask = cell(1, K);
    for j = 1:K, mask{j} = true(1, mdds.domain(j)); end
    mask{l} = false(1, d);
    mask{l}(k) = true;
    mass(k) = mdd_rec_masked(mdds, g, mask);
end
end
