function g = ssa_nrm_guards(sn)
% G = SSA_NRM_GUARDS(SN)
% The per-feature tests the NRM engine's dispatch and gate both read.
%
% ONE BODY, THREE CALLERS. SOLVER_SSA_ANALYZER reads the fields INDIVIDUALLY,
% because its explicit 'nrm' arm falls back to the serial engine with a warning
% that names the offending construct, so it has to know WHICH test failed.
% SSA_NRM_ELIGIBLE reads them together to answer "should the NRM be preferred".
% SSA_NRM_REFUSAL reads only `sched`, because that is the one condition the
% explicit arm does not fall back on.
%
% These were local functions of SOLVER_SSA_ANALYZER until 2026-09-05, which is
% why the eligibility test could not be asked from a gate. Moving them here is
% what lets all three callers share one definition instead of three copies.
%
% Fields (all true = the NRM can take it):
%   sched   every station's discipline has a reaction form
%   phase   non-exponential service only where every job present is in service
%   renege  no station renegs with non-exponential patience
%   balk    no balking strategy outside QUEUE_LENGTH
%   gd      no global (Whittle) dependence
%   cache   every cache configuration (always true; see below)
%   fcr     every finite capacity region rule (always true; see below)
%   spn     a Petri net the NRM's reaction builder can read (see below)

% see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch, EXT-source exclusion)
allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
    SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS, ...
    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
    SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING];
g.sched = all(arrayfun(@(s) any(s == allowedSched), sn.sched));

% Phase expansion splits the class-level share across a class's phases in the
% ratio kir/nir, which needs only the per-phase populations -- true of the
% INF/PS family, where every job present is in service. A buffered policy
% instead needs the phase multiset of the jobs ACTUALLY in service, which the
% waiting-only buffer does not record, so non-exponential service there still
% needs the serial engine.
exact = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.LPS, ...
    SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, ...
    SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT];
g.phase = true;
for ist = 1:sn.nstations
    for r = 1:sn.nclasses
        if sn.procid(ist,r) == ProcessType.DISABLED || sn.procid(ist,r) == ProcessType.EXP
            continue
        end
        if ~any(sn.sched(ist) == exact)
            g.phase = false;
            break
        end
    end
    if ~g.phase
        break
    end
end

% The NRM abandons at the aggregate rate (waiting count)*mu, which is only
% correct when patience is memoryless; phase-type patience would need each
% waiting job's remaining phase. SOLVER_SSA rejects the same combination.
g.renege = true;
if isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass)
    bad = (sn.impatienceClass == ImpatienceType.RENEGING) & (sn.impatienceType ~= ProcessType.EXP);
    g.renege = ~any(bad(:));
end

% QUEUE_LENGTH is a pure function of the state vector, so the NRM draws it at
% firing time; EXPECTED_WAIT and COMBINED depend on the mean waiting time and
% need the serial engine (State.afterEventStation rejects them likewise).
g.balk = true;
if isfield(sn,'balkingStrategy') && ~isempty(sn.balkingStrategy)
    bs = sn.balkingStrategy(:);
    g.balk = all(bs == 0 | bs == BalkingStrategy.QUEUE_LENGTH);
end

% The NRM builds one propensity closure per reaction from the per-station
% population slice; a global handle reads the whole population matrix, which
% that closure does not receive. SOLVER_SSA (serial) carries the factor.
g.gd = ~isfield(sn,'gdscaling') || isempty(sn.gdscaling);

% True for every Cache node. The NRM models a cache access as an immediate
% state-dependent class switch (read -> hit/miss/retrieval) at the cache node,
% applying the same replacement logic as State.afterEventCache to the cache
% contents carried alongside the buffers, INCLUDING the retrieval (delayed-hit)
% system: a miss for an item not yet being fetched begins a retrieval and a
% concurrent request for an item already being fetched is absorbed as a delayed
% hit -- matching the serial engine's sample-path semantics.
g.cache = true;

% Finite capacity regions are supported under both rules. DROP destroys a
% refused job; WAITQ parks it in a per-region FIFO and admits it head-of-line as
% capacity frees. The NRM carries the FIFO explicitly (see fcrReleaseCascade),
% so no region rule forces a fallback. Linear-constraint and memory-budget
% regions ride the same admission test.
g.fcr = true;

% The NRM SPN path (solver_ssa_nrm, spnBuildMode) builds one constant-rate
% reaction per timed mode and per Source arrival, so it has no form for a
% non-exponential timed firing, a non-exponential Source arrival into a Place,
% a Source that reaches no Place, or an infinite initial marking. Each of these
% it RAISED from inside the builder, which the default dispatch reached
% unconditionally on every net; the serial engine (State.afterGlobalEvent)
% serves all four, so they are a fallback condition and not a solver refusal.
% A marking-dependent firing rate is NOT among them: no SSA engine applies the
% g(marking) multiplier, and SSA_FIRINGDEP_REFUSAL refuses it under every
% method before an engine is chosen.
g.spn = true;
if any(sn.nodetype == NodeType.Transition)
    R = sn.nclasses;
    for ind = 1:sn.nnodes
        switch sn.nodetype(ind)
            case NodeType.Transition
                np = sn.nodeparam{ind};
                for m = 1:np.nmodes
                    if np.timing(m) ~= TimingStrategy.IMMEDIATE
                        fK = np.firingphases(m);
                        if isnan(fK) || fK ~= 1 || isempty(np.firingproc{m})
                            g.spn = false;
                            return
                        end
                    end
                end
            case NodeType.Source
                ist = sn.nodeToStation(ind);
                for r = 1:R
                    lambda = sn.rates(ist, r);
                    if isnan(lambda) || lambda <= 0
                        continue
                    end
                    if sn.procid(ist, r) ~= ProcessType.EXP
                        g.spn = false;
                        return
                    end
                    feeds = false;
                    for jnd = find(sn.nodetype == NodeType.Place)'
                        if any(sn.rtnodes((ind-1)*R + r, (jnd-1)*R + (1:R)) > 0)
                            feeds = true;
                            break
                        end
                    end
                    if ~feeds
                        g.spn = false;
                        return
                    end
                end
            case NodeType.Place
                isf = sn.nodeToStateful(ind);
                if sn.isstateful(ind) && isfield(sn,'state') && numel(sn.state) >= isf ...
                        && ~isempty(sn.state{isf})
                    [~, nir] = State.toMarginalAggr(sn, ind, sn.state{isf});
                    if any(isinf(nir))
                        g.spn = false;
                        return
                    end
                end
        end
    end
end
end
