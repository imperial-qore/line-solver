function [bool, reason] = fluid_qsys_admits(sn, method)
% [BOOL, REASON] = FLUID_QSYS_ADMITS(SN, METHOD)
%
% @brief The shape a single-station fluid limit is stated for, asked as a gate.
%
% 'ggisgi.fluid', 'ggingi.tga', 'tvms', 'mtginf' and 'mol' are closed forms for
% ONE open class through ONE Source and ONE queueing station, not integrations
% of the network drift, so they cannot answer any other shape. Three of them
% are limits for a queue customers ABANDON and need a reneging patience law;
% three need a finite server count (the Erlang formula of 'mol', the sqrt(n)
% fluctuation of 'tga', the staffing of 'tvms'). Every rule here used to be an
% inline error inside SOLVER_FLUID_QSYS_ANALYZER, invisible to a caller: a
% report offered 'mol' on a Source -> Delay -> Sink model and the run stopped on
% the infinite server count.
%
% A time-varying SERVICE law is refused too: the analyzer reads the service ccdf
% off a stationary MAP pair and a schedule slot (NHPP/MAPt/PHt) is not one, so
% it silently fell back to an exponential of the same mean. The ARRIVAL side is
% where a schedule belongs, and SN_ARRIVAL_RATE_FUN reads it there.
%
% The registry cannot name any of this ("exactly one station", "requires a
% patience law", "a finite server count"), hence the structural predicate.
% Called by SOLVER_FLUID_QSYS_ANALYZER, so the run stops on it, and by
% FLUID_METHOD_REFUSAL, so a report sees the same verdict. One predicate, two
% callers. The horizon rule of the three time-varying limits is a separate
% predicate, FLUID_QSYS_HORIZON, because it is a rule on the OPTIONS.
%
% @param sn NetworkStruct of the model
% @param method the method name, any spelling SolverFLD.canonicalMethod takes
% @return bool true when METHOD may run on this model
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
reason = '';
m = SolverFLD.canonicalMethod(method);
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queueing = sn.nodeToStation(sn.nodetype == NodeType.Queue | sn.nodetype == NodeType.Delay);
if sn.nclasses ~= 1 || numel(source_ist) ~= 1 || numel(queueing) ~= 1 || sn.nclosedjobs > 0
    reason = sprintf(['The ''%s'' method is a single-station limit: it needs one open class ' ...
        'through one Source and one queueing station.'], m);
    return
end
queue_ist = queueing(1);
if any(strcmp(m, {'ggisgi.fluid','ggingi.tga','tvms'}))
    if isempty(sn_patience_handles(sn, queue_ist, 1))
        reason = sprintf(['The ''%s'' method needs a reneging patience law on the queue ' ...
            '(Queue.setPatience): it is a limit for a queue customers abandon.'], m);
        return
    end
end
if any(strcmp(m, {'ggingi.tga','tvms','mol'}))
    s = sn.nservers(queue_ist);
    if ~isfinite(s) || s < 1
        reason = sprintf('The ''%s'' method needs a finite number of servers.', m);
        return
    end
end
if isfield(sn, 'procid') && ~isempty(sn.procid)
    pid = sn.procid(queue_ist, 1);
    if pid == ProcessType.NHPP || pid == ProcessType.MAPT || pid == ProcessType.PHT
        reason = sprintf(['The ''%s'' method is stated for a stationary service law, and the ' ...
            'queue''s service process is a time-varying schedule; put the schedule on the ' ...
            'arrival process instead.'], m);
        return
    end
end
bool = true;
end
