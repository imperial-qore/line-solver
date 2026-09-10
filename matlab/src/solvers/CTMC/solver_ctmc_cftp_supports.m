function [bool, reason] = solver_ctmc_cftp_supports(sn, options)
% [BOOL, REASON] = SOLVER_CTMC_CFTP_SUPPORTS(SN, OPTIONS)
%
% @brief Can the cftp perfect sampler be asked for this model?
%
% The model-class gate of the cftp method, asked as a predicate rather than
% raised. SOLVER_CTMC_CFTP calls it to refuse a model it cannot sample, and
% SolverCTMC.supportsModelMethod calls it so that a CALLER (model.help,
% findSolver, SolverAUTO) sees the same verdict before paying for a run.
% Keeping it in one place is the point: a second copy in the analyzer is how
% the report and the run drift into two different answers.
%
% The sampler is exact only on the closed single-class product form its
% balance function encodes; anything else must be refused, not approximated.
% What the feature registry CAN name (open classes, non-product-form
% disciplines, non-exponential service, ...) is also declared in
% SolverCTMC.getMethodFeatureSet. The structural rules the registry has no name
% for are the class count, the station count and the steady-state restriction.
% The rest is stated TWICE ON PURPOSE: SolverCTMC.supportsModelMethod asks this
% predicate BEFORE the feature gate, and the analyzer asks it INSTEAD of the
% feature gate, so a construct declared only in the featset would be sampled
% away by a run with options.enableChecks=false.
%
% @param sn NetworkStruct of the model
% @param options solver options (read for timespan)
% @return bool true when the cftp sampler may run
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
if ~isinf(options.timespan(1))
    reason = 'The cftp method supports steady-state analysis only, not transient analysis.';
    return
end
if sn.nclasses ~= 1
    reason = sprintf('The cftp method supports single-class models only, this model has %d classes.',sn.nclasses);
    return
end
if any(isinf(sn.njobs)) || sn.njobs(1) < 1
    reason = 'The cftp method supports closed models only, with a finite positive population.';
    return
end
if sn.nstations < 2
    reason = 'The cftp method requires at least two stations.';
    return
end
for ind=1:sn.nnodes
    switch sn.nodetype(ind)
        case {NodeType.Queue, NodeType.Delay, NodeType.Router}
            % product-form stations and stateless routers only
        otherwise
            reason = sprintf('The cftp method supports Queue, Delay and Router nodes only, node %d is of a different type.',ind);
            return
    end
end
pfSched = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.FCFS, ...
    SchedStrategy.SIRO, SchedStrategy.LCFSPR];
for i=1:sn.nstations
    if ~any(sn.sched(i) == pfSched)
        reason = sprintf('The cftp method requires a product-form scheduling strategy (INF, PS, FCFS, SIRO, LCFSPR) at station %d.',i);
        return
    end
    if sn.phases(i,1) > 1
        reason = sprintf('The cftp method requires exponential service times, station %d has %d phases.',i,sn.phases(i,1));
        return
    end
    if isfinite(sn.cap(i)) && sn.cap(i) < sn.njobs(1)
        reason = sprintf('The cftp method requires infinite buffers, station %d has capacity %d.',i,sn.cap(i));
        return
    end
    if ~isfinite(sn.rates(i,1)) || sn.rates(i,1) <= 0
        reason = sprintf('The cftp method requires a finite positive service rate at station %d.',i);
        return
    end
end
% The four State constructs the Kijima-Matsui balance function has no term
% for. Each is withdrawn by name in SolverCTMC.getMethodFeatureSet's cftp
% branch and repeated here because SOLVER_CTMC_CFTP asks this predicate and not
% the featset: without them the sampler would draw from the plain Gordon-Newell
% chain and report the answer as if the construct were absent.
if isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass) && any(sn.impatienceClass(:) ~= 0)
    reason = 'The cftp method does not support impatience (reneging or balking): the balance function it samples has no abandonment.';
    return
end
if isfield(sn,'balkingStrategy') && ~isempty(sn.balkingStrategy) && any(sn.balkingStrategy(:) ~= 0)
    reason = 'The cftp method does not support balking: the balance function it samples admits every arrival.';
    return
end
if isfield(sn,'retrialProc') && ~isempty(sn.retrialProc) && any(~cellfun(@isempty, sn.retrialProc(:)))
    reason = 'The cftp method does not support retrial orbits: the balance function it samples has no orbit.';
    return
end
if isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) && any(sn.hasbreakdown(:) == 1)
    reason = 'The cftp method does not support server breakdowns: the balance function it samples has no outage.';
    return
end
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    reason = 'The cftp method does not support load-dependent, class-dependent or joint-dependent service rates.';
    return
end
if isfield(sn,'nregions') && sn.nregions > 0
    reason = 'The cftp method does not support finite capacity regions.';
    return
end
pfRouting = [RoutingStrategy.PROB, RoutingStrategy.RAND, RoutingStrategy.DISABLED];
for ind=1:sn.nnodes
    for r=1:sn.nclasses
        if ~any(sn.routing(ind,r) == pfRouting)
            reason = sprintf('The cftp method requires Markovian routing (PROB, RAND), node %d uses a state-dependent strategy.',ind);
            return
        end
    end
end
bool = true;
reason = '';
end
