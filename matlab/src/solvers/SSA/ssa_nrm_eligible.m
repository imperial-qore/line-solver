function [ok, reason] = ssa_nrm_eligible(sn)
% [OK, REASON] = SSA_NRM_ELIGIBLE(SN)
% Whether the next-reaction-method engine can run this model AT ALL.
%
% THIS IS THE PREFERENCE QUESTION, not the gate question. SOLVER_SSA_ANALYZER
% asks it to decide whether the 'default' and 'parallel' paths should run the
% NRM instead of the serial engine. A caller asking whether an EXPLICIT
% options.method='nrm' will produce an answer must ask SSA_NRM_REFUSAL instead,
% which is narrower: the explicit arm FALLS BACK to the serial engine, with a
% warning, for every one of these conditions except the scheduling one.
%
% The tests themselves live in SSA_NRM_GUARDS, which SOLVER_SSA_ANALYZER also
% reads field by field so its fallback warnings can name the construct. One
% body, three callers.
%
% The NRM simulates open and closed models alike -- it only lacks Fork/Join node
% handling -- so the model-class exclusion is Fork/Join, not the INF/PS-only
% sn_is_population_model.

g = ssa_nrm_guards(sn);
ok = g.sched && g.phase && g.renege && g.balk && g.gd && g.cache && g.fcr ...
    && g.spn && ~sn_has_fork_join(sn);
reason = '';
if ~ok
    reason = ['The ''nrm'' engine cannot run this model; the serial engine can. ' ...
        'Use options.method=''serial''.'];
end
end
