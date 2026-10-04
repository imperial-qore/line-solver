function G = mdd_rec(mdds, g)
% G = MDD_REC(MDDS, G_L)
% Normalising constant of a product-form model whose reachable set is held in a
% decision diagram: G = sum_{s in S} prod_l g_l(s_l), by Algorithm 1 of the
% MDD-rec paper.
%
% FORMALISM-AGNOSTIC. Nothing here knows what a level is: on the lattice
% sum_k s_k = n of a closed queueing network this collapses to Buzen's
% convolution (paper Appendix B), and on an S-invariant reachable Petri net to
% the Coleman-Henderson-Taylor convolution, SPN_CONV (Sec. 5). Unlike either, it
% needs only that the reachable set be finite and encoded -- no lattice, no
% S-invariant reachability.
%
% -- Input
% MDDS : MDD.toStruct of the reachable set
% G_L  : 1 x K cell, G_L{l}(v+1) is g_l(v)
% -- Output
% G    : the normalising constant
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 4.
%
% See also MDD_REC_MASKED, MDD_REC_MARGINAL, SPN_METRICS, SPN_CONV.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

G = mdd_rec_masked(mdds, g, []);
end
