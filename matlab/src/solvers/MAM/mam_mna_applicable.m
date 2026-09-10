function [tf, reason] = mam_mna_applicable(sn)
% [TF, REASON] = MAM_MNA_APPLICABLE(SN)
%
% Can the 'mna' method of SolverMAM answer this model? REASON is '' when it
% can and otherwise names what does not fit.
%
% These are the CORRECTNESS rules of SOLVER_MNA_OPEN and SOLVER_MNA_CLOSED,
% the ones a by-name request must clear; the default-method chooser in
% solver_mam_analyzer adds its own PREFERENCES on top (a multiserver PS
% station, a multiclass FCFS station), which are about accuracy and cost and
% do not refuse a caller who asks for mna by name.
%
%   - no mixed open/closed model (the two analyzers are one each);
%   - round-robin routing in open models only (the deterministic split is
%     carried by the open traffic equations; the closed sweep has none);
%   - no self-looping closed class (no inter-station flow to decompose);
%   - no fork-join (both analyzers raise 'Fork nodes not supported yet');
%   - a closed model with one class per chain: SOLVER_MNA_CLOSED bisects over
%     CLASSES but stores the throughput in the CHAIN-indexed lambda and
%     renormalizes chain c with the class-indexed sn.njobs(c), which is only
%     the same quantity when the two indexings coincide;
%   - INF, PS, FCFS and EXT stations only: a station under any other
%     discipline is never updated by the flow sweep and keeps a zero queue
%     length, which the table then reports as the answer.
%
% ONE PREDICATE, THREE CALLERS: SolverMAM.supportsModelMethod, the 'mna' arm
% of solver_mam_analyzer, and its default-method chooser (mnaApplies).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
isOpen = sn_is_open_model(sn);
isClosed = sn_is_closed_model(sn);
if ~isOpen && ~isClosed
    reason = 'The mna method does not support mixed open/closed models.';
    return
end
if ~isOpen && any(sn.routing(:) == RoutingStrategy.RROBIN)
    reason = 'The mna method supports round-robin routing in open models only.';
    return
end
if sn_has_fork_join(sn)
    reason = 'The mna method does not support fork-join (Fork nodes not supported yet by the QNA sweep). Use the dec.source method.';
    return
end
if isClosed && sn.nchains ~= sn.nclasses
    reason = sprintf(['The mna method supports a closed model only when every chain holds one ' ...
        'class (this model has %d classes in %d chains): solver_mna_closed indexes the ' ...
        'throughput by chain and the population by class. Use the dec.source method.'], ...
        sn.nclasses, sn.nchains);
    return
end
for ist = 1:sn.nstations
    switch sn.sched(ist)
        case {SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.EXT}
            % covered by the flow sweep
        otherwise
            reason = sprintf(['The mna method does not support the %s scheduling strategy at ' ...
                'station %d: the flow sweep updates INF, PS and FCFS stations only and would ' ...
                'report a zero queue length there. Use the dec.source method.'], ...
                SchedStrategy.toText(sn.sched(ist)), ist);
            return
    end
end
% A self-looping closed class (confined to a single non-INF station) has no
% inter-station flow to decompose; see _kb/06-solver-catalog.md for rationale
V = cellsum(sn.visits);
for k = 1:sn.nclasses
    if ~isfinite(sn.njobs(k))
        continue;
    end
    vis = find(V(:, k) > GlobalConstants.FineTol);
    if isscalar(vis) && sn.sched(vis) ~= SchedStrategy.INF && sn.sched(vis) ~= SchedStrategy.EXT
        reason = sprintf(['The mna method does not support self-looping ' ...
            'classes (class %d is confined to station %d with no ' ...
            'inter-station flow to decompose). Use the dec.source method.'], k, vis);
        return
    end
end
tf = true;
reason = '';
end
