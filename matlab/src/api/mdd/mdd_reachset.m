function mdd = mdd_reachset(domain, init, nextfun)
% MDD = MDD_REACHSET(DOMAIN, INIT, NEXTFUN)
% Generate and store the reachability set into a quasi-reduced ordered MDD,
% using the decision diagram itself as the visited-state store, after
% A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and Storage
% Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.
%
% -- Input
% DOMAIN  : 1 x K vector of per-level local-state counts (values 0..domain(k)-1)
% INIT    : 1 x K initial global state (0-based local values)
% NEXTFUN : function handle s -> T, returning an (m x K) matrix whose rows are
%           the successor states of s under the next-state function
% -- Output
% MDD     : an MDD object holding every state reachable from INIT
%
% -- Remarks
% This is the basic (explicit-frontier) realisation: a breadth-first search
% enumerates successors while the MDD provides the O(K) membership test that
% replaces the usual explicit visited hash. The stored set lives entirely in
% the MDD (O(number of nodes) memory); only the transient BFS frontier is held
% explicitly. Symbolic image computation / saturation, which removes the
% explicit frontier too, is the natural next step but is out of scope here.
%
% See also: MDD, mdd_closedqn.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

mdd = MDD(domain);
mdd.insert(init);

frontier = init;          % rows are states still to be expanded
head = 1;
while head <= size(frontier, 1)
    s = frontier(head, :);
    head = head + 1;
    T = nextfun(s);
    for r = 1:size(T, 1)
        t = T(r, :);
        if ~mdd.member(t)
            mdd.insert(t);
            frontier(end + 1, :) = t; %#ok<AGROW>
        end
    end
    % drop already-expanded rows periodically to bound frontier memory
    if head > 1024 && 2 * head > size(frontier, 1)
        frontier = frontier(head:end, :);
        head = 1;
    end
end
mdd.compact();   % reclaim the dead nodes left by the append-only inserts
end
