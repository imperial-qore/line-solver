function [ok, reason] = lqns_method_refusal(model, method)
% [OK, REASON] = LQNS_METHOD_REFUSAL(MODEL, METHOD)
% Whether lqns or lqsim can serve this LayeredNetwork under METHOD, as a
% predicate. METHOD '' asks the method-neutral half only.
%
% ONE PREDICATE, TWO CALLERS. SolverLQNS.supportsModelMethod asks it, so
% model.help and SolverAUTO never offer a pair that dies at run time, and
% runAnalyzer asks it again before the .lqnx is written, so a caller naming the
% method by hand gets LINE's own sentence rather than the binary's parse error.
%
% WHY NOT A FEATURE SET. SolverLQNS.supports used to compare the per-layer
% Networks that SolverLN builds against a flat set naming Queue, Exp and FCFS.
% Every LN layer carries a Delay('Clients') under INF with PROB routing, which
% that set did not declare, so the comparison refused every layered model there
% is, and it never saw an LQN-level construct at all. lqns reads the .lqnx and
% not the layers, so the rules are stated on the LayeredNetwork itself.
%
% THE RULES. (1) The constructs neither binary models, named one by one by
% SolverLQNS.unsupportedLNConstructs. (2) The scheduling vocabulary: writeXML
% writes SchedStrategy.toText verbatim (FCFSPRPRIO as 'pri'), and the LQN schema
% names fcfs, ps, inf, hol, pri and ref only, so any other discipline reaches
% lqns as a syntax error for what is a modelling limit. (3) Under 'sim'/'lqsim',
% a replicated processor or task, which lqsim refuses to simulate.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = true;
reason = '';
if nargin < 2 || isempty(method)
    method = '';
end
method = char(method);

reasons = SolverLQNS.unsupportedLNConstructs(model);
if ~isempty(reasons)
    ok = false;
    reason = sprintf(['SolverLQNS cannot serve this model: %s. Use SolverLN ' ...
        '(the native layered solver), which models all three.'], strjoin(reasons, '; '));
    return
end

lqnxSched = {'fcfs','ps','inf','hol','pri','ref'};
elems = [model.hosts(:); model.tasks(:)];
for k = 1:numel(elems)
    elem = elems{k};
    if ~isa(elem, 'LayeredNetworkElement') || ~isprop(elem, 'scheduling') || isempty(elem.scheduling)
        continue
    end
    schedText = char(elem.scheduling);
    if strcmpi(schedText, SchedStrategy.toText(SchedStrategy.FCFSPRPRIO))
        schedText = 'pri'; % the spelling writeXML emits for it
    end
    if ~any(strcmpi(schedText, lqnxSched))
        ok = false;
        reason = sprintf(['''%s'' is scheduled %s, which the LQN XML schema does not name: lqns ' ...
            'and lqsim accept fcfs, ps, inf, hol, pri (preemptive priority) and ref only. Use ' ...
            'SolverLN or SolverLDES.'], elem.name, upper(schedText));
        return
    end
end

if any(strcmpi(method, {'sim','lqsim'}))
    for k = 1:numel(elems)
        elem = elems{k};
        if isa(elem, 'LayeredNetworkElement') && isprop(elem, 'replication') ...
                && ~isempty(elem.replication) && elem.replication > 1
            ok = false;
            reason = sprintf(['''%s'' is replicated %d times, which lqsim refuses to simulate. ' ...
                'Use an lqns method, SolverLN or SolverLDES, which materialise the replicas.'], ...
                elem.name, elem.replication);
            return
        end
    end
end
end
