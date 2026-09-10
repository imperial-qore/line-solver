function [bool, reason, isAoI] = fluid_mfq_admits(sn)
% [BOOL, REASON, ISAOI] = FLUID_MFQ_ADMITS(SN)
%
% @brief Can the 'mfq' method (and its aliases 'butools', 'aoi') run here?
%
% The Markovian fluid queue answers ONE open station, Source -> Queue -> Sink.
% SOLVER_FLUID_ANALYZER used to fall back to the matrix method on any other
% shape, and SOLVER_MFQ_PRIO did the same inside its own branch, so a caller's
% 'mfq' label was answered by a different algorithm under that label. A listed
% name must run as itself, so the shapes are refused by name instead.
%
% Three arms, tried in the analyzer's own order:
%   the age-of-information shape (AOI_IS_AOI): one open class, one server, a
%     buffer of 1 or 2, FCFS/LCFS/LCFSPR. The finite buffer IS the model there,
%     so ISAOI lets a caller exempt it from the binding-capacity gate, as 'mol'
%     is exempt for the same reason.
%   the priority branch (SOLVER_MFQ_PRIO), taken when the classes carry
%     distinct priorities: a class-independent service rate, a MAP {D0,D1}
%     arrival for every open class, at least one of them modulated. Each of
%     those used to be a silent fallback inside the branch.
%   the plain branch (SOLVER_MFQ), which reads the arrival and service
%     processes of class 1 only, so it is stated for ONE open class.
% Closed classes are refused by the feature set (SolverFLD.getMethodFeatureSet)
% rather than here: SOLVER_MFQ skips them and reports zeros.
%
% Called by SOLVER_FLUID_ANALYZER, so the run stops on it, and by
% FLUID_METHOD_REFUSAL, so a report sees the same verdict before running. One
% predicate, two callers.
%
% @param sn NetworkStruct of the model
% @return bool true when 'mfq' may run on this model
% @return reason the refusal, or '' when BOOL is true
% @return isAoI true when the age-of-information arm is the one that runs

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
reason = '';
isAoI = aoi_is_aoi(sn);
if isAoI
    bool = true;
    return
end
[ok, info] = fluid_is_single_queue(sn);
if ~ok
    reason = sprintf(['The ''mfq'' method is a single-queue fluid model (Source -> Queue -> Sink, ' ...
        'one or infinitely many servers): %s.'], info.errorMsg);
    return
end
openClasses = find(isinf(sn.njobs));
% the same test SOLVER_FLUID_ANALYZER switches on
if numel(unique(sn.classprio)) > 1
    qi = info.queueStation;
    si = info.sourceStation;
    mu = sn.rates(qi, openClasses);
    if isempty(openClasses) || any(~isfinite(mu)) || any(mu <= 0)
        reason = ['The priority branch of ''mfq'' needs a finite positive service rate for ' ...
            'every open class.'];
        return
    end
    if any(abs(mu - mu(1)) > GlobalConstants.FineTol * max(1, mu(1)))
        reason = ['The priority branch of ''mfq'' is a fluid priority queue drained at one ' ...
            'rate, so it needs a class-independent service rate.'];
        return
    end
    nph = 1;
    for k = openClasses
        proc = sn.proc{si}{k};
        if ~iscell(proc) || numel(proc) < 2
            reason = sprintf(['The priority branch of ''mfq'' needs a MAP {D0,D1} arrival process ' ...
                'for class %d.'], k);
            return
        end
        nph = nph * size(proc{1}, 1);
    end
    if nph < 2
        reason = ['The priority branch of ''mfq'' needs at least one Markov-modulated (multi-phase) ' ...
            'arrival process: with exponential arrivals its fluid model degenerates.'];
        return
    end
    bool = true;
    return
end
if numel(openClasses) ~= 1
    reason = sprintf(['The ''mfq'' method reads the arrival and service processes of one open ' ...
        'class and this model has %d; several classes are served only by its priority branch ' ...
        '(distinct class priorities).'], numel(openClasses));
    return
end
bool = true;
end
