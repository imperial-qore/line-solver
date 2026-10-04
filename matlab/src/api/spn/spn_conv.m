function G = spn_conv(S, V, g)
% G = SPN_CONV(S, V, G_L)
% G = SPN_CONV(INV, G_L)
% Normalising constant of an S-invariant reachable product-form stochastic
% Petri net, by convolution over the invariant load vector.
%
% With S the minimal-support S-invariant matrix and V = S m0 the load vector,
% the reachability set of an S-INVARIANT REACHABLE net is exactly
% {m >= 0 : S m = V}, and conditioning on the marking of one place partitions it
% (Lemma 5.1 of the MDD-rec paper). Writing G_j(W) for the mass of the markings
% supported on the first j places with S m = W,
%
%   G_0(W) = [W == 0],   G_j(W) = sum_i g_j(i) G_{j-1}(W - i S_j),
%
% and G = G_n(V). On a net whose only invariant is "the tokens are conserved"
% this is Buzen's convolution for a closed queueing network, one place per
% station.
%
% NO ILP IS SOLVED. Coleman-Henderson-Taylor obtain the marking set M_p(P',W)
% from the feasibility of an integer program so that the sum skips the terms
% that contribute nothing. Here the sum simply runs over i whose residual
% W - i S_j stays non-negative and the recursion returns zero on an infeasible
% residual, which gives the same value: the ILP is an optimisation of the
% enumeration, not part of the definition. Memoising on (j, W) keeps the walk
% over the reachable residuals rather than over all of them.
%
% S-INVARIANT REACHABILITY IS NOT CHECKED, and cannot be cheaply: no algorithm
% is known that decides it without generating the reachability set. On a net
% that fails it, {m : S m = V} is strictly larger than the reachable set and
% this returns a normalising constant over unreachable markings too, which is
% why MDD_REC -- which walks the reachable set itself -- is the general
% algorithm and this one the special case. Compare the two on a new net before
% trusting this one on it.
%
% -- Input
% S   : S(i,p), the minimal-support S-invariants, one row per invariant; or the
%       struct SPN_SINVARIANTS returned, in which case V holds the factors
% V   : the load vector S m0, one entry per invariant
% G_L : 1 x n cell, G_L{p}(i+1) is g_p(i), the product-form factor of i tokens
%       in place level p; its length bounds the marking of that level
% -- Output
% G   : the normalising constant
%
% -- Reference
% J. Coleman, W. Henderson, P. Taylor, "Product form equilibrium distributions
% and a convolution algorithm for stochastic Petri nets", Performance
% Evaluation 26(3), 1996, 159-180; the point of comparison for MDD-rec in
% S. Balsamo, A. Marin, I. Stojic, Future Generation Computer Systems 111
% (2020) 475-490, Sec. 5.1.
%
% See also MDD_REC, SPN_SINVARIANTS, SPN_MDD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isstruct(S)                      % spn_conv(inv, g) off the invariant basis
    g = V;
    V = S.V;
    S = S.S;
end
if isempty(S)
    line_error(mfilename, 'the net has no S-invariant to convolve over');
end
if size(S, 1) ~= numel(V)
    line_error(mfilename, 'one load-vector entry per invariant is required');
end
n = numel(g);
if size(S, 2) ~= n
    line_error(mfilename, 'the invariant matrix and g must agree on the place-level count');
end
if any(S(:) < 0)
    line_error(mfilename, ['an S-invariant has a negative weight, so the residual recursion ' ...
        'has no monotone bound on the marking']);
end

memo = dictionary(string.empty, []);
G = i_rec(n, V(:)');

    function val = i_rec(j, W)
        if j == 0
            val = double(all(W == 0));
            return
        end
        key = sprintf('%d|%s', j, sprintf('%d,', W));
        if isKey(memo, key)
            val = memo(key);
            return
        end
        p = j;
        acc = 0;
        for i = 0:numel(g{p}) - 1
            rem = W - i * S(:, p)';
            if any(rem < 0)
                break                % S is non-negative, so larger i only gets worse
            end
            if g{p}(i + 1) == 0, continue; end
            acc = acc + g{p}(i + 1) * i_rec(j - 1, rem);
        end
        memo(key) = acc;
        val = acc;
    end
end
