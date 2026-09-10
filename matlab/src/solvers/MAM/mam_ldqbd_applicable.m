function [tf, reason, regime] = mam_ldqbd_applicable(sn)
% [TF, REASON, REGIME] = MAM_LDQBD_APPLICABLE(SN)
%
% Can the 'ldqbd' method of SolverMAM answer this model? REASON is '' when it
% can and otherwise names what does not fit; REGIME is 'closed' (one Delay
% and one FCFS Queue, finite population), 'open' (one Source and one FCFS
% Queue, Poisson arrivals) or '' on a refusal.
%
% THE RULE IS A SHAPE, which a feature set cannot state: SOLVER_MAM_LDQBD
% builds one level-dependent QBD whose level is the queue length of THE
% station, so it needs exactly two stations and a single class, and the
% open regime needs a Poisson stream because a MAP would need an arrival
% phase the chain does not carry. The stability test of the open regime is
% left to the analyzer: it depends on the load-dependent capacity factor the
% analyzer derives, not on the model shape.
%
% ONE PREDICATE, FOUR CALLERS. SolverMAM.supportsModelMethod asks it for
% 'ldqbd', and for 'default' on a load-dependent model (the closed regime is
% the only MAM path that reads sn.lldscaling); solver_mam_analyzer asks it to
% route a single-class closed Delay+Queue to the exact chain; solver_mam_ldqbd
% asks it before building anything. Until 2026-09-05 each carried its own
% copy of the shape and the report offered ldqbd on models the analyzer then
% refused by name.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
regime = '';
if sn.nclasses ~= 1
    reason = 'The ldqbd method requires a single-class model.';
    return
end
nDelay  = sum(sn.sched == SchedStrategy.INF);
nQueue  = sum(sn.sched == SchedStrategy.FCFS);
nSource = sum(sn.sched == SchedStrategy.EXT);
if ~isfinite(sn.njobs(1))
    if nSource ~= 1 || nQueue ~= 1 || sn.nstations ~= 2
        reason = 'Open LDQBD method requires exactly one Source and one Queue station.';
        return
    end
    srcIdx = find(sn.sched == SchedStrategy.EXT, 1);
    % External Poisson arrivals only (MAP/MMPP arrivals are not yet supported).
    if sn.procid(srcIdx, 1) ~= ProcessType.EXP || (isfield(sn,'phases') && ~isempty(sn.phases) ...
            && sn.phases(srcIdx, 1) > 1)
        reason = ['Open LDQBD method currently supports Poisson (exponential) ' ...
            'arrivals only; the Source uses a MAP/MMPP process.'];
        return
    end
    regime = 'open';
else
    if sn.njobs(1) <= 0
        reason = 'Closed LDQBD method requires a positive population.';
        return
    end
    if nDelay ~= 1 || nQueue ~= 1 || sn.nstations ~= 2
        reason = 'Closed LDQBD method requires exactly one Delay and one Queue station.';
        return
    end
    regime = 'closed';
end
tf = true;
reason = '';
end
