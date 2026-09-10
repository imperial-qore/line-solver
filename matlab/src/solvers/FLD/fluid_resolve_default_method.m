function [method, reason] = fluid_resolve_default_method(sn, options, model)
% [METHOD, REASON] = FLUID_RESOLVE_DEFAULT_METHOD(SN, OPTIONS, MODEL)
%
% Which concrete method options.method='default' stands for, for lang='matlab'.
%
% Preference order: 'rmf' for cache models, then 'minnormal' whenever
% FLUID_MINNORMAL_APPLICABLE accepts the model, then the historical choice of
% 'closing' for DPS and 'matrix' otherwise. The second-order closure dominates
% the first-order methods on every family measured against exact CTMC and is
% the only method that can represent GPS at all, so it is preferred wherever it
% applies.
%
% This must be called BEFORE the feature gate, not from the dispatch switch:
% NetworkSolver.runAnalyzerChecks validates getMethodFeatureSet(options.method),
% and 'default' is not 'minnormal', so a GPS model would be rejected by the
% gate before ever reaching the resolution. Keeping the decision in one
% function is what stops the gate and the dispatch from disagreeing, which
% would silently send a GPS model to the matrix method.
%
% Parameters:
%   sn      - NetworkStruct, after sn_nonmarkov_toph
%   options - solver options
%   model   - the Network, optional. Only the binding-capacity test needs it,
%             because it reads the caps the USER set rather than the derived
%             sn.classcap (see NetworkSolver.checkBindingCapacity); omitting it
%             skips that branch, leaving the caller with the pre-existing
%             preference order.
%
% Returns:
%   method - the resolved method name
%   reason - why 'minnormal' was declined, empty when it was selected
%
% See also FLUID_MINNORMAL_APPLICABLE, SOLVER_FLUID_MOMENTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
% A stochastic Petri net has exactly one fluid route, SOLVER_FLUID_PETRI, and it
% is reached through 'dae': the marking closure is stated as a differential-
% algebraic system because a P-invariant and an immediate firing flow are both
% equations rather than events. No other method has a drift for a Place at all.
if any(sn.nodetype == NodeType.Transition)
    method = 'dae';
    reason = 'the model is a stochastic Petri net';
    return
end
if any(sn.nodetype == NodeType.Cache)
    method = 'rmf';
    reason = 'the model has cache nodes';
    return
end
% A BINDING BUFFER OR A CAPACITY REGION ALSO HAS ONE FLUID ROUTE, for the same
% reason a Petri net does: nothing else in the FLD tree reads sn.cap, sn.classcap
% or the region limit, so every other method integrates the capped station as an
% unbounded one -- which is exactly why RUNANALYZER refuses them. Resolving
% 'default' to one of those turned a model this solver CAN answer into an error
% whose advice was to type the very method the resolution should have picked
% ('SolverFLD(model)' on cqn_bas_blocking refused, 'SolverFLD(model,''method'',
% ''dae'')' solved it). The test is the gate's own, so the two cannot disagree.
%
% Only when the dae route actually accepts the model: where it does not, falling
% through leaves the refusal to the gate, which names the blocking feature.
if nargin >= 3 && isa(model, 'Network')
    hasRegion = isfield(sn,'nregions') && ~isempty(sn.nregions) && sn.nregions > 0;
    capBinds = ~NetworkSolver.checkBindingCapacity(model, 'SolverFLD');
    if hasRegion || capBinds
        daeOk = fluid_dae_applicable(sn, options);
        if daeOk
            method = 'dae';
            if hasRegion
                reason = 'the model has a finite capacity region';
            else
                reason = 'the model has a binding finite buffer';
            end
            return
        end
        % Why dae declined is deliberately not carried further: 'minnormal' is
        % asked next and the gate below then refuses the model naming the
        % capacity itself, which is the blocker the caller has to act on.
    end
end
[ok, reason] = fluid_minnormal_applicable(sn, options);
if ok
    method = 'minnormal';
    return
end
if sn_has_dps(sn)
    method = 'closing';
else
    method = 'matrix';
end
end
