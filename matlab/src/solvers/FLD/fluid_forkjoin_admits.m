function [bool, reason] = fluid_forkjoin_admits(sn, method)
% [BOOL, REASON] = FLUID_FORKJOIN_ADMITS(SN, METHOD)
%
% @brief Can METHOD run the fluid fork-join fixed point on this model?
%
% A fork-join model is not integrated as one drift: the MMT transform replaces
% the fork by auxiliary OPEN classes arriving at a Source it adds (closed model
% or not, see ModelAdapter.mmt) and the answer is the fixed point of solving
% that mixed model repeatedly (@NetworkSolver/fjFixedPoint). So the inner solve
% is always a mixed network, and a method refused on one is refused here by
% name rather than by a failure inside the fixed point: 'softmin' and
% 'statedep' have no EXT branch, 'refined' is closed-only, and on an OPEN model
% the DAE form has no unknowns for the auxiliary classes, so its inner solve
% fails on the class count rather than returning a drift.
%
% Called by @SolverFLD/runAnalyzer, so the run stops on it, and by
% @SolverFLD/supportsModelMethod, so a CALLER sees the same verdict before
% paying for the fixed point. One predicate, two callers.
%
% @param sn NetworkStruct of the model
% @param method the concrete method name
% @return bool true when METHOD may run the fork-join fixed point here
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = true;
reason = '';
if ~any(sn.nodetype == NodeType.Fork)
    return
end
m = SolverFLD.canonicalMethod(method);
switch m
    case {'softmin','statedep'}
        % THE TRANSFORM ALWAYS ADDS A SOURCE, closed model or not: ModelAdapter.mmt
        % creates one for the auxiliary open classes that carry the parallelism
        % (rate GlobalConstants.FineTol). ODE_SOFTMIN and ODE_STATEDEP have no
        % branch for an EXT station and raise on it, so the inner solve of every
        % fork-join model died there while the closed-model featset admitted it.
        bool = false;
        reason = sprintf(['The %s method has no route through the fork-join fixed point: the MMT ' ...
            'transform hands the inner solve a mixed network with a Source, and its ODE has no ' ...
            'branch for one (it refuses open models). Use options.method=''closing''.'], m);
    case 'refined'
        % The inner solve is that same mixed network, on which 'refined' is
        % refused by its own closed-model rule (see SolverFLD.getMethodFeatureSet):
        % the 1/N correction is solved over the full state, source pool included.
        bool = false;
        reason = ['The refined method has no route through the fork-join fixed point: the MMT ' ...
            'transform hands the inner solve a mixed network, and refined is a closed-model ' ...
            'method. Use options.method=''minnormal'', the same closure without the 1/N term.'];
    case 'dae'
        if ~any(isinf(sn.njobs))
            return
        end
        bool = false;
        reason = ['The dae method has no route through the fork-join fixed point on an OPEN model: ' ...
            'the MMT transform hands the inner solve a mixed network whose auxiliary open classes ' ...
            'the DAE form carries no unknowns for. Use options.method=''minnormal'', which is the ' ...
            'same closure and does run that fixed point.'];
end
end
