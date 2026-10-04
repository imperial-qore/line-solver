function [ok, reason] = fluid_minnormal_applicable(sn, options)
% [OK, REASON] = FLUID_MINNORMAL_APPLICABLE(SN, OPTIONS)
%
% Whether SOLVER_FLUID_MOMENTS can answer this model, used by the 'default'
% method of SolverFLD to prefer 'minnormal' over 'matrix' when it applies.
%
% The test is STATIC: it inspects the model and the options, never the
% solution, so a feature the model declares is decided here and only here.
% The one condition that cannot be static is a NON-HYPERBOLIC fluid fixed
% point (balanced bottlenecks, a saturated multiclass station, an overloaded
% open station): it exists only once the mean is solved. FLUID_LYAPUNOV
% detects it and raises 'LINE:FluidNonHyperbolic', and
% @SolverFLD/runAnalyzer switches a RESOLVED 'minnormal' to the first-order
% method on that identifier alone. An explicit options.method='minnormal'
% still fails loudly, so the closure never absorbs a real defect silently.
%
% Every condition below mirrors a guard that SOLVER_FLUID_MOMENTS,
% FLUID_MOMENT_TERMS or @SolverFLD/runAnalyzer would otherwise raise, so this
% function and those guards must move together.
%
% Parameters:
%   sn      - NetworkStruct, after sn_nonmarkov_toph
%   options - solver options
%
% Returns:
%   ok     - true when 'minnormal' can be selected
%   reason - one line naming the blocking feature, empty when OK
%
% See also SOLVER_FLUID_MOMENTS, FLUID_MOMENT_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = false;

% Open and mixed models are supported: FLUID_MOMENT_TERMS projects the EXT
% source pool out of the covariance. What it cannot take is a NON-POISSON
% arrival stream, whose source coordinates track the phase of a single arrival
% process rather than a population.
for ist = 1:sn.nstations
    if sn.sched(ist) ~= SchedStrategy.EXT
        continue
    end
    for r = 1:sn.nclasses
        if ~isnan(sn.rates(ist,r)) && ~isempty(sn.mu{ist}{r}) && numel(sn.mu{ist}{r}) > 1
            reason = sprintf('class %d has a %d-phase (non-Poisson) arrival process', ...
                r, numel(sn.mu{ist}{r}));
            return
        end
    end
end
% A cache model is answered through the decomposition analyzer, with the
% closure in its network step, so cache nodes no longer decline the method.
% What still declines is a replacement strategy with no drift-based fluid
% model, mirroring the runtime guard in SOLVER_FLD_CACHEQN_ANALYZER: LRU,
% HLRU, CLIMB and QLRU are answered by a characteristic-time fixed point,
% which is not a fluid method and has no covariance.
caches = find(sn.nodetype == NodeType.Cache);
for ci = 1:numel(caches)
    ch = sn.nodeparam{caches(ci)};
    if ~isfield(ch,'replacestrat') || ~(ch.replacestrat == ReplacementStrategy.RR ...
            || ch.replacestrat == ReplacementStrategy.FIFO ...
            || ch.replacestrat == ReplacementStrategy.SFIFO)
        reason = 'a cache uses a replacement strategy with no drift-based fluid model';
        return
    end
end

% scheduling policies with a branch in ODE_RATES_CLOSING_FACTORS; anything
% else falls through to rates = x, i.e. an infinite server
supported = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
    SchedStrategy.FCFS, SchedStrategy.DPS, SchedStrategy.GPS];
for ist = 1:sn.nstations
    if ~any(sn.sched(ist) == supported)
        reason = sprintf('station %d uses %s, which has no fluid drift branch', ...
            ist, SchedStrategy.toText(sn.sched(ist)));
        return
    end
end

% the Lyapunov solve is cubic in the phase-resolved state, so the same cap
% SOLVER_FLUID_MOMENTS enforces decides selection rather than being hit later
maxstate = 200;
if isfield(options,'config') && isfield(options.config,'moment_maxstate') ...
        && ~isempty(options.config.moment_maxstate)
    maxstate = options.config.moment_maxstate;
end
nstate = 0;
for ist = 1:sn.nstations
    for r = 1:sn.nclasses
        if ~isnan(sn.rates(ist,r)) && ~isempty(sn.mu{ist}{r}) && ~any(isnan(sn.mu{ist}{r}))
            nstate = nstate + numel(sn.mu{ist}{r});
        end
    end
end
if nstate > maxstate
    reason = sprintf('the phase-resolved state has %d coordinates, above the moment_maxstate limit of %d', ...
        nstate, maxstate);
    return
end

% the moment methods need the untransformed event set and an autonomous drift
if isfield(options,'config') && isfield(options.config,'rate_traj') && ~isempty(options.config.rate_traj)
    reason = 'options.config.rate_traj makes the drift time-varying';
    return
end
if isfield(options,'config') && isfield(options.config,'nhpp_sched') && ~isempty(options.config.nhpp_sched)
    reason = 'options.config.nhpp_sched makes the drift time-varying';
    return
end

ok = true;
reason = '';
end
