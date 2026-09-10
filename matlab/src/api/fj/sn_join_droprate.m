function DropRateJoin = sn_join_droprate(sn, TN, AN)
% DROPRATEJOIN = SN_JOIN_DROPRATE(SN, TN, AN)
%
% Rate at which sibling tasks are discarded at each Join, as an
% (nstations x nclasses) matrix that is zero away from the Join rows.
%
% A Join is the one station where the loss identity ArvR - Tput does NOT hold,
% because the two rates are in different units: AN counts the SIBLINGS offered
% to the join (N per parent job) while TN counts the PARENT jobs released by
% it (one per synchronisation). Reading ArvR - Tput there reports (N-1)/N of
% the offered traffic as lost at every join, standard joins included, when a
% standard join loses nothing at all.
%
% The siblings a join actually consumes are K per synchronisation, where K is
% the quorum (K = N on a standard join), so
%
%   DropRateJoin = max(0, AN - K*TN)
%
% which is 0 for a standard join and (N-K)*TN for a quorum, the rate at which
% the stragglers of an already-fired parent are discarded on arrival.
%
% This is the DERIVED value, exact given TN and AN. A solver that MEASURES the
% discards on its own sample path (SolverLDES) reports its own; see
% NetworkSolver.getAvgLossTable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

DropRateJoin = zeros(sn.nstations, sn.nclasses);
if isempty(TN) || isempty(AN) || ~isfield(sn,'fj') || isempty(sn.fj) || ~any(sn.fj(:))
    return
end
for j = find(sn.nodetype == NodeType.Join)'
    ist = sn.nodeToStation(j);
    if ist < 1 || ist > size(DropRateJoin,1)
        continue
    end
    for r = 1:sn.nclasses
        a = AN(ist,r);
        t = TN(ist,r);
        if ~isfinite(a) || ~isfinite(t) || a <= 0
            continue
        end
        % PER CLASS: a variable forking level makes the sibling count differ
        % between classes, so it cannot be hoisted out of this loop.
        nsib = sn_join_siblings(sn, j, r);
        if nsib <= 0
            continue
        end
        kreq = sn_join_quorum(sn, j, r, nsib);
        DropRateJoin(ist,r) = max(0, a - kreq*t);
    end
end
end
