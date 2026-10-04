function [ok, reason] = ssa_firingdep_refusal(sn)
% [OK, REASON] = SSA_FIRINGDEP_REFUSAL(SN)
% Whether SolverSSA can answer a Petri net whose firing rates depend on the
% marking (Transition.setFiringRateDependence).
%
% IT CANNOT, UNDER ANY METHOD, AND SAYS SO BY NAME. The NRM builds one
% constant-propensity reaction per timed mode, so the g(marking) multiplier
% has no place in it, and the contract the solver is held to (line-test
% test_spn_firing_dependence, test_ssa_rejects_dependence) is a refusal
% rather than a nominal-rate answer. The registry has no name for the
% construct ('Transition' and 'Firing' are declared and the rule is about
% the handle behind them), so it is structural.
%
% ONE PREDICATE, THREE CALLERS: SolverSSA.supportsModelMethod answers with
% it so model.help does not offer any ssa.* row on such a net,
% solver_ssa_analyzer raises with it before an engine is chosen, and the
% NRM reaction builder raises with it on the enableChecks=false path.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = true;
reason = '';
for ind = 1:sn.nnodes
    if sn.nodetype(ind) ~= NodeType.Transition
        continue
    end
    np = sn.nodeparam{ind};
    if ~isfield(np, 'firingdep') || isempty(np.firingdep)
        continue
    end
    for m = 1:min(np.nmodes, numel(np.firingdep))
        if ~isempty(np.firingdep{m})
            ok = false;
            reason = sprintf(['Transition %s mode %d uses a marking-dependent firing rate ' ...
                '(setFiringRateDependence), which SolverSSA does not support; use SolverCTMC ' ...
                'or SolverLDES.'], sn.nodenames{ind}, m);
            return
        end
    end
end
end
