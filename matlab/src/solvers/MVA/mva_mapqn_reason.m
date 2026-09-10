function reason = mva_mapqn_reason(sn)
% REASON = MVA_MAPQN_REASON(SN)
%
% The structural premise of the SolverMVA method 'amva.mapqn': the reason it
% cannot solve SN, or '' when it can. One predicate, three callers:
% listValidMethods drops the name when it is nonempty, supportsModelMethod
% reports it, and solver_mva_mapqn_analyzer raises it, so the offered,
% reported and run answers cannot drift apart.
%
% The premise is the model of mapqn_amva: a closed network of exactly one
% infinite-server station with exponential think times and one FCFS
% single-server queue whose service processes are Markovian (any process
% with a (D0,D1) representation: Exp, Erlang, HyperExp, PH, APH, Coxian,
% MAP, MMPP2), every class cycling delay -> queue -> delay. Only the
% station count, the server count, the think-time law and the routing shape
% are structural; MAP support itself is a feature-set matter and lives in
% SolverMVA.getMethodFeatureSet.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
tol = 1e-12;
if sn.nstations ~= 2
    reason = 'Method ''amva.mapqn'' requires exactly two stations: one delay (infinite server) and one FCFS queue.';
    return
end
isDelay = (sn.sched(:) == SchedStrategy.INF) | isinf(sn.nservers(:));
if sum(isDelay) ~= 1
    reason = 'Method ''amva.mapqn'' requires exactly one delay (infinite-server) station and one queue.';
    return
end
id = find(isDelay); iq = find(~isDelay);
if sn.sched(iq) ~= SchedStrategy.FCFS
    reason = sprintf('Method ''amva.mapqn'' requires FCFS scheduling at the queue; station %d is %s.', iq, SchedStrategy.toText(sn.sched(iq)));
    return
end
if sn.nservers(iq) ~= 1
    reason = 'Method ''amva.mapqn'' supports a single-server queue only.';
    return
end
if any(isinf(sn.njobs)) || sn.nclosedjobs <= 0
    reason = 'Method ''amva.mapqn'' supports closed models only.';
    return
end
markovian = [ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.PH, ...
    ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ProcessType.MMPP2];
R = sn.nclasses;
sd = sn.stationToStateful(id); sq = sn.stationToStateful(iq);
for r = 1:R
    if sn.njobs(r) <= 0
        continue
    end
    if sn.procid(id, r) ~= ProcessType.EXP
        reason = sprintf('Method ''amva.mapqn'' requires exponential think times; class %d has a %s think time.', r, ProcessType.toText(sn.procid(id, r)));
        return
    end
    if ~any(sn.procid(iq, r) == markovian) || isempty(sn.proc{iq}{r})
        reason = sprintf('Method ''amva.mapqn'' requires a Markovian (MAP-representable) service process at the queue; class %d is not.', r);
        return
    end
    if abs(sn.rt((sd - 1) * R + r, (sq - 1) * R + r) - 1) > tol || abs(sn.rt((sq - 1) * R + r, (sd - 1) * R + r) - 1) > tol
        reason = sprintf('Method ''amva.mapqn'' requires every class to cycle delay -> queue -> delay without class switching; class %d does not.', r);
        return
    end
end
end
