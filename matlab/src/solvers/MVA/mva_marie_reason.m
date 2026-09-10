function reason = mva_marie_reason(sn)
% REASON = MVA_MARIE_REASON(SN)
%
% The structural premise of the SolverMVA method 'marie': the reason it cannot
% solve SN, or '' when it can. One predicate, two callers: supportsModelMethod
% reports it (through SolverMVA.supportsMarie) and solver_mva_marie_analyzer
% raises it, so the reported and the run answers cannot drift apart.
% listValidMethods withholds the name on the two halves that make the whole
% family pointless, an open model and class-dependent routing; the other
% halves are reported here rather than hidden.
%
% The premise is the model of pfqn_marie: a closed network whose
% infinite-server stations fold into a per-chain think time and whose queueing
% stations are FCFS (service-sensitive, the SCV is read), PS or LCFS-PR
% (insensitive, taken exponential), each isolated from its chain demands alone,
% so no cache, no fork-join and no load-, class- or joint-dependent scaling,
% with a multiserver station admitted only when there is a single chain, the
% one case the isolation carries a server count. The nameable halves (open
% class, cache, fork-join, scaling, discipline) are also dropped from the
% method's feature set in SolverMVA.getMethodFeatureSet, so the report names
% the offending feature; the rest is structural and lives here alone.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if any(isinf(sn.njobs)) || any(sn.nodetype == NodeType.Source)
    reason = ['The ''marie'' method supports closed models only; this model has open ' ...
        'classes. Use another SolverMVA method (e.g. ''default'').'];
    return
end
if sn_has_classdep_routing(sn)
    reason = ['The ''marie'' method aggregates the classes of a chain into one demand ' ...
        'vector, which is exact only when every class is routed alike; this model ' ...
        'switches class or routes its classes with different probabilities. Use ' ...
        'another SolverMVA method (e.g. ''default'').'];
    return
end
if any(sn.nodetype == NodeType.Cache) || sn_has_fork_join(sn)
    reason = ['The ''marie'' method has no cache or fork-join arm. Use another ' ...
        'SolverMVA method (e.g. ''default'').'];
    return
end
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    reason = ['The ''marie'' method isolates each station from its chain demands and ' ...
        'SCVs alone and carries no load-, class- or joint-dependent scaling. Use ' ...
        'another SolverMVA method (e.g. ''default'').'];
    return
end
isDelay = isinf(sn.nservers(:)) | (sn.sched(:) == SchedStrategy.INF);
schedOK = isDelay | (sn.sched(:) == SchedStrategy.FCFS) | ...
    (sn.sched(:) == SchedStrategy.PS) | (sn.sched(:) == SchedStrategy.LCFSPR);
if ~all(schedOK)
    bad = find(~schedOK, 1);
    reason = sprintf(['The ''marie'' method supports FCFS, PS, LCFSPR and Delay ' ...
        'stations only; station %d has an unsupported scheduling strategy. Use ' ...
        'another SolverMVA method.'], bad);
    return
end
nservers = sn.nservers(~isDelay);
nservers(~isfinite(nservers)) = 1;
if sn.nchains > 1 && any(nservers > 1)
    reason = ['The ''marie'' method supports multiserver queueing stations for ' ...
        'single-chain models only; this model is multichain with a multiserver station.'];
end
end
