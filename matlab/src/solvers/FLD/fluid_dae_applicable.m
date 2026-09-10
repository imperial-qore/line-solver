function [ok, reason] = fluid_dae_applicable(sn, options)
% [OK, REASON] = FLUID_DAE_APPLICABLE(SN, OPTIONS)
%
% Whether SOLVER_FLUID_DAE can answer this model, used by @SolverFLD/runAnalyzer
% to try 'dae' before dropping a declined 'minnormal' to a first-order method.
%
% WHY THIS EXISTS SEPARATELY FROM FLUID_MINNORMAL_APPLICABLE. The two methods
% state the SAME closure and differ only in how the coupled equations are
% discharged, so a model 'minnormal' accepts is almost always one 'dae' accepts
% too. Almost: 'dae' carries a finite-difference Jacobian over the whole
% unknown vector rather than one Lyapunov solve, so its state cap is lower; it
% closes on the per-station variance only, so DPS and GPS are out; and it has
% no decomposition route, so a cache model is out. Those three are exactly the
% difference set, and naming them here keeps the fallback ladder from entering
% a rung that would refuse the model a moment later.
%
% The test is STATIC, for the same reason FLUID_MINNORMAL_APPLICABLE is: a rung
% chosen by trial and rollback would make the reported method depend on a failed
% run. The one condition that cannot be static is the NON-HYPERBOLIC fixed point
% the ladder exists to route around -- it exists only once the mean is solved --
% and 'dae' fails on it loudly, which is what moves the ladder to its last rung.
%
% Every condition below mirrors a refusal SOLVER_FLUID_DAE would otherwise
% raise, so this function and those refusals must move together.
%
% Parameters:
%   sn      - NetworkStruct, after sn_nonmarkov_toph
%   options - solver options
%
% Returns:
%   ok     - true when 'dae' can be selected
%   reason - one line naming the blocking feature, empty when OK
%
% See also SOLVER_FLUID_DAE, FLUID_MINNORMAL_APPLICABLE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = false;

% The same closure as 'minnormal', so the same arrival-stream rule: FLUID_MOMENT_-
% TERMS raises on a multi-phase source, whose coordinates track the phase of one
% arrival process rather than a population and carry no linear noise
% approximation. Without this line the ladder in runAnalyzer entered the dae
% rung on such a model and the rung refused it a moment later.
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

% A cache model is a DECOMPOSITION, and 'dae' has no arm for it: SOLVER_FLUID_-
% ANALYZER routes 'minnormal' on a cache model through SOLVER_FLD_CACHEQN_-
% ANALYZER, with the closure inside the network step, and sends 'dae' straight
% to SOLVER_FLUID_DAE, which would read the cache nodes as ordinary stations.
if any(sn.nodetype == NodeType.Cache)
    reason = 'a cache model is answered by the decomposition analyzer, which has no dae route';
    return
end

% The DPS and GPS shares close on the covariance BETWEEN station coordinates,
% not on the station total, so their closure state is a matrix block rather than
% the scalar the Newton vector carries. Mirrors the refusal in SOLVER_FLUID_DAE.
for ist = 1:sn.nstations
    if sn.sched(ist) == SchedStrategy.DPS || sn.sched(ist) == SchedStrategy.GPS
        reason = sprintf(['station %d uses %s, whose share closes on the covariance between its ' ...
            'class coordinates rather than on the station variance'], ...
            ist, SchedStrategy.toText(sn.sched(ist)));
        return
    end
end

% The simultaneous solve is quartic overall, against the cubic of one Lyapunov
% solve, so its crossover is lower than the 200 SOLVER_FLUID_MOMENTS permits and
% it carries its own cap. Same count as FLUID_MINNORMAL_APPLICABLE forms, so the
% two limits are read on the same scale.
maxstate = 100;
if isfield(options,'config') && isfield(options.config,'dae_maxstate') ...
        && ~isempty(options.config.dae_maxstate)
    maxstate = options.config.dae_maxstate;
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
    reason = sprintf(['the phase-resolved state has %d coordinates, above the dae_maxstate limit ' ...
        'of %d'], nstate, maxstate);
    return
end

ok = true;
reason = '';
end
