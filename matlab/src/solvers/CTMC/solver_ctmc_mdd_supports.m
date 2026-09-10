function [bool, reason, P, kind] = solver_ctmc_mdd_supports(sn, options)
% [BOOL, REASON, P, KIND] = SOLVER_CTMC_MDD_SUPPORTS(SN, OPTIONS)
%
% @brief Can the mdd decision-diagram method be asked for this model?
%
% The model-shape gate of the mdd method, asked as a predicate rather than
% raised. SOLVER_CTMC_MDD_ANALYZER calls it before it builds anything, and
% SolverCTMC.supportsModelMethod calls it so that a CALLER (model.help,
% findSolver, SolverAUTO) sees the same verdict without paying for a run.
% One predicate with two callers is what stops the report and the analyzer
% from disagreeing about which models the method serves.
%
% A STOCHASTIC PETRI NET IS EXEMPT from the shape rules: a Place model is
% read through SPN_MDD, which builds the reachable set and the Kronecker
% descriptor from the marking rather than from the (station,class) encoding,
% so neither the single-class rule nor the closed-population rule applies to
% it. The transient rule binds it all the same.
%
% THE TWO DEEPER RULES ARE HERE TOO, and the analyzer reads their products
% off this predicate rather than recomputing them: P is the station-to-station
% routing chain of the single class, which must be stochastic (a completion
% must move the job to another station, so a Router or a leak refuses), and
% KIND is the local-state encoding the disciplines and service laws admit,
% 'np' (count, or count plus one non-preemptive phase) or 'ps' (per-phase
% counts, shared servers only). Both are decided from sn alone, so a caller
% is told about them at the same price as the class count. On a net, or on a
% refusal, P and KIND are empty.
%
% @param sn NetworkStruct of the model
% @param options solver options (read for timespan); optional
% @return bool true when the mdd method may run
% @return reason the refusal, or '' when BOOL is true
% @return P station-to-station routing matrix of the single class
% @return kind 'np' or 'ps', the local-state encoding to build

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
reason = '';
P = [];
kind = '';
% The level iteration answers a stationary question only: it never forms
% the generator uniformization would need, so a finite horizon has nothing
% to integrate. RUNANALYZER used to route a transient request here and
% return the steady state under it.
if nargin >= 2 && ~isempty(options) && isfield(options,'timespan') ...
        && ~isempty(options.timespan) && ~isinf(options.timespan(1))
    reason = 'the mdd method supports steady-state analysis only, not transient analysis';
    return
end
if any(sn.nodetype == NodeType.Place)
    bool = true;
    return
end
if sn.nclasses ~= 1
    reason = sprintf(['the mdd method analyses single-class networks; this model ' ...
        'has %d classes. The Kronecker descriptor would need one level per (station,class).'], sn.nclasses);
    return
end
if any(sn.nodetype == NodeType.Source) || any(sn.nodetype == NodeType.Sink)
    reason = ['the mdd method analyses CLOSED networks; an open stream makes the ' ...
        'marking unbounded, so the reachable set has no finite decision diagram'];
    return
end
if ~isfinite(sn.njobs(1)) || sn.njobs(1) <= 0
    reason = 'the mdd method needs a finite positive closed population';
    return
end

% ---- the station-to-station routing chain of the single class, from sn.rt
% (which is indexed over stateful nodes, class-major)
M = sn.nstations;
R = sn.nclasses;
P = zeros(M);
for i = 1:M
    ni = sn.stationToStateful(i);
    for j = 1:M
        nj = sn.stationToStateful(j);
        P(i, j) = sn.rt((ni - 1) * R + 1, (nj - 1) * R + 1);
    end
end
rs = sum(P, 2);
if any(abs(rs - 1) > 1e-8)
    P = [];
    reason = ['the station-to-station routing chain is not stochastic; the mdd ' ...
        'method needs every completion to move the job to another station'];
    return
end

% ---- which local-state encoding represents these disciplines exactly.
% Exponential service is discipline-insensitive for the queue-length law, so
% the compact count encoding serves any work-conserving station. Phase-type
% service is not: the count-plus-one-phase encoding is non-preemptive, while
% processor sharing needs the per-phase counts of every job present.
nph = ones(1, M);
if isfield(sn, 'phases') && ~isempty(sn.phases)
    nph = sn.phases(:, 1)';
end
phAt = find(nph > 1);
if isempty(phAt)
    kind = 'np';
    bool = true;
    return
end
isShared = false(1, M); isNP = false(1, M);
for i = 1:M
    s = sn.sched(i);
    isShared(i) = ismember(s, [SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
        SchedStrategy.INF]);
    isNP(i) = ismember(s, [SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
        SchedStrategy.HOL]) && sn.nservers(i) == 1;
end
if all(isShared)
    kind = 'ps';
    bool = true;
    return
end
if all(isNP(phAt))
    kind = 'np';
    bool = true;
    return
end
bad = phAt(~isNP(phAt) & ~isShared(phAt));
if isempty(bad), bad = phAt(1); end
P = [];
reason = sprintf(['station %d combines a phase-type service law with a discipline ' ...
    'that neither local encoding represents: the count-plus-phase encoding is non-preemptive, ' ...
    'and the per-phase-count encoding covers only shared servers (PS/DPS/GPS/INF). Mixing a ' ...
    'shared and a non-preemptive phase-type station in one model is likewise unsupported.'], bad(1));
end
