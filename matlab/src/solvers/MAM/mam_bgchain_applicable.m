function [tf, reason] = mam_bgchain_applicable(sn, options)
% [TF, REASON] = MAM_BGCHAIN_APPLICABLE(SN, OPTIONS)
%
% Can the 'bgchain' method of SolverMAM answer this model? REASON is '' when
% it can and otherwise names what is missing.
%
% The rules are the ones SOLVER_MAM_BGCHAIN and MAM_BGCHAIN_CTMC raise:
%   - at least one closed chain (the background chain IS the closed
%     population vector, so a purely open model has nothing to build it from);
%   - the closed chains visit some station (the support of the chain);
%   - no class priorities at an HOL or FCFSPRPRIO station (the open classes
%     are aggregated into one phase-type mixture per station, which cannot
%     express a priority order);
%   - no fork-join (the chain conserves the closed population per station,
%     which a fork violates);
%   - a background chain within options.config.bgstates_max states, sized by
%     MAM_BGCHAIN_STATES without building it.
%
% ONE PREDICATE, THREE CALLERS. SolverMAM.supportsModelMethod asks it so the
% method is not offered on a model it cannot answer; solver_mam_analyzer asks
% it before choosing bgchain as the closed or mixed default; solver_mam_bgchain
% asks it before building anything, so a caller naming the method gets the
% same sentence. Until 2026-09-05 the analyzer kept a private copy
% (bgchainApplies) and the gate a second one, and the state-space cap sat in
% neither, so the report offered bgchain on models MAM_BGCHAIN_CTMC refused.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = struct();
end
tf = false;

if sn_is_open_model(sn)
    reason = ['The bgchain method requires at least one closed class: the ' ...
        'background chain IS the closed population vector, which a purely open ' ...
        'model does not have. Use the dec.source method.'];
    return
end
if any(sn.sched == SchedStrategy.HOL | sn.sched == SchedStrategy.FCFSPRPRIO) ...
        && any(sn.classprio ~= sn.classprio(1))
    reason = ['The bgchain method does not support class priorities: it ' ...
        'aggregates the open classes into one phase-type mixture per station, ' ...
        'which cannot express a priority order. Use the dec.source method.'];
    return
end
if sn_has_fork_join(sn)
    reason = ['The bgchain method does not support fork-join: the background ' ...
        'chain conserves the closed population per station, which a fork ' ...
        'violates. Use the dec.source method.'];
    return
end
[~, ~, Vchain, ~, Nchain] = sn_get_demands_chain(sn);
closedChains = find(isfinite(Nchain) & Nchain > 0);
if isempty(closedChains)
    reason = ['The bgchain method requires at least one closed class: the ' ...
        'background chain IS the closed population vector, which a purely open ' ...
        'model does not have. Use the dec.source method.'];
    return
end
if ~any(any(Vchain(:, closedChains) > GlobalConstants.Zero))
    reason = 'The closed classes of this model visit no station.';
    return
end
% The background chain enumerates the closed population vector over the
% stations the closed classes visit; the cap is the one MAM_BGCHAIN_CTMC
% enforces, so the gate refuses exactly what the run would.
if isfield(options, 'config') && isfield(options.config, 'bgstates_max') ...
        && ~isempty(options.config.bgstates_max)
    bgstates_max = options.config.bgstates_max;
else
    bgstates_max = 20000;
end
nstates = mam_bgchain_states(sn, options);
if nstates > bgstates_max
    reason = sprintf(['The background chain of this model has %d states, above the ' ...
        'options.config.bgstates_max = %d cap of the bgchain method. Raise the cap or ' ...
        'use the dec.source method.'], nstates, bgstates_max);
    return
end
tf = true;
reason = '';
end
