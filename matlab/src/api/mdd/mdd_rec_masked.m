function G = mdd_rec_masked(mdds, g, mask)
% G = MDD_REC_MASKED(MDDS, G_L)
% G = MDD_REC_MASKED(MDDS, G_L, MASK)
% Unnormalised mass of the masked subset of a reachable set held in a decision
% diagram, i.e. Algorithm 1 of the MDD-rec paper under a per-level restriction.
%
% A product-form model has P(s) = (1/G) prod_k g_k(s_k) over its levels, and
%
%   G = sum_{s in S} prod_k g_k(s_k).
%
% Summing state by state is exponential and numerically unstable. MDD-rec
% instead walks the diagram that already encodes S, accumulating the
% unnormalised mass of each node ONCE (Def. 4.4, Algorithm 1):
%
%   M(<l.p>) = sum_{v in S_l} g_l(v) * M(<l.p>[v]),  M(TRUE)=1, M(FALSE)=0
%
% so the cost is O(sum_l |nodes_l| * |S_l|) rather than O(|S|), and
% G = M(root).
%
% THE MASK is how Sec. 5.3 of the paper computes measures. Restricting the sum
% at level l to a subset of its local values gives the unnormalised mass of the
% corresponding subset of S, so P(m_l = k) and P(e_j >= k) are the same
% recursion under a different mask rather than three separate algorithms.
%
% -- Input
% MDDS : MDD.toStruct of the reachable set (level 1 is the root)
% G_L  : 1 x K cell, G_L{l}(v+1) is g_l(v), the per-level product-form factor
% MASK : 1 x K cell of logical row vectors; [] admits everything, which is the
%        plain MDD-rec of Algorithm 1 and returns the normalising constant
% -- Output
% G    : the unnormalised mass of the masked subset
%
% -- Note
% The g_l themselves, and the test that the model has a product form at all,
% are the caller's: the paper declares that out of scope (Sec. 3.2). Passing
% g_l that do not describe a product-form model returns a number that is not
% the normalising constant of anything, and nothing here can detect it.
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 4.
%
% See also MDD_REC, MDD_REC_MARGINAL, SPN_REC_ENABLED, MDD_REACHSET.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3, mask = []; end
K = mdds.K;
if numel(g) ~= K
    line_error(mfilename, 'one g_l per level is required');
end
for l = 1:K
    if numel(g{l}) ~= mdds.domain(l)
        line_error(mfilename, 'g_l must have one entry per local state');
    end
end
if ~isempty(mask)
    if numel(mask) ~= K
        line_error(mfilename, 'the mask must have one row per level');
    end
    for l = 1:K
        if numel(mask{l}) ~= mdds.domain(l)
            line_error(mfilename, 'the mask must have one entry per local state');
        end
    end
end

if mdds.root == MDD.TERM_FALSE
    G = 0;
    return
end

memo = cell(1, K);
for l = 1:K
    memo{l} = zeros(mdds.nnodes(l), 1);
end

% Levels are visited bottom-up, so every child mass a node needs is already
% memoised when the node is reached; an explicit sweep avoids the recursion
% depth a level count in the hundreds would otherwise reach.
for l = K:-1:1
    for id = 1:mdds.nnodes(l)
        arcs = mdds.node{l}(id, :);
        acc = 0;
        for v = 1:mdds.domain(l)
            if ~isempty(mask) && ~mask{l}(v), continue; end
            gv = g{l}(v);
            if gv == 0, continue; end
            ch = arcs(v);
            if l == K
                if ch ~= MDD.TERM_TRUE, continue; end
                acc = acc + gv;
            else
                if ch == MDD.TERM_FALSE, continue; end
                acc = acc + gv * memo{l + 1}(ch);
            end
        end
        memo{l}(id) = acc;
    end
end
G = memo{1}(mdds.root);
end
