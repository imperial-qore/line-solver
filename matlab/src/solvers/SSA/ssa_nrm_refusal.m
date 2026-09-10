function [ok, reason] = ssa_nrm_refusal(sn)
% [OK, REASON] = SSA_NRM_REFUSAL(SN)
% Whether an EXPLICIT options.method='nrm' can answer this model.
%
% NOT SSA_NRM_ELIGIBLE, and the difference is the whole point of having two
% predicates. SSA_NRM_ELIGIBLE asks whether the NRM should be PREFERRED, and it
% is deliberately wide: SOLVER_SSA_ANALYZER consults it on the 'default' and
% 'parallel' paths to pick an engine. This one asks what a GATE must ask, which
% is whether the name the caller typed will produce an answer.
%
% The two differ because the explicit 'nrm' arm of SOLVER_SSA_ANALYZER FALLS
% BACK to the serial engine, with a warning, for reneging patience, balking,
% phase-type service, global dependence and cache retrieval. A model carrying
% any of those still gets a correct answer under 'nrm', so refusing it here
% would withdraw a row the solver honours -- the report would be wrong in the
% opposite direction.
%
% What is left is the ONE condition that reaches SOLVER_SSA_ANALYZER_NRM and
% raises there rather than falling back: a scheduling discipline the reaction
% network has no form for (solver_ssa_analyzer_nrm.m:24,
% 'solver_ssa_analyzer_nrm:UnsupportedPolicy'). That is the rule this gate
% states, and only that one.
%
% The C++ port gates on the wider set ON PURPOSE and is not diverging: its NRM
% has no fallback arm at all, so every one of the six conditions raises there.
% Each gate states its own engine's reach.

ok = true;
reason = '';

% The one field of SSA_NRM_GUARDS the explicit arm does not fall back on.
g = ssa_nrm_guards(sn);
if ~g.sched
    % see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch, EXT-source exclusion)
    allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
        SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS, ...
        SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
        SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
        SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
        SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING];
    unsupported = unique(sn.sched(~arrayfun(@(s) any(s == allowedSched), sn.sched)));
    names = strjoin(arrayfun(@(s) SchedStrategy.toText(s), unsupported, ...
        'UniformOutput', false), ', ');
    ok = false;
    reason = sprintf(['The ''nrm'' engine has no reaction form for the %s scheduling ' ...
        'policy. Use options.method=''serial''.'], names);
end
end
