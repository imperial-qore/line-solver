function [ok, reason] = ldes_ln_refusal(model)
% [OK, REASON] = LDES_LN_REFUSAL(MODEL)
% Whether the LDES layered engine serves this LayeredNetwork, as a predicate.
%
% The engine is the JVM backend (jline.solvers.ldes), which validates the LQN
% against its SolverLDES.getLNFeatureSet at run time; the MATLAB
% SolverLDES.getLNFeatureSet is the mirror of that declaration. supports() used
% to answer true for every LayeredNetwork and leave the refusal to Java, so
% model.help offered 'ldes' on an LQN with a DPS processor or a Zipf demand and
% the run then died in the backend. Two callers: SolverLDES.supports (hence
% supportsModelMethod, model.help and SolverAUTO) and runAnalyzer before the
% Java run, so both speak the same sentence.
%
% WHAT IS COMPARED. The LayeredNetwork has no feature recorder of its own (its
% getUsedLangFeatures returns the per-layer Networks SolverLN builds, which
% carry no LQN-level name), so the two dimensions that vary across LQNs are
% marked here: the discipline of every processor and task (SchedStrategy_*) and
% the law of every host demand and think time (the distribution's feature name,
% the one the Network recorder marks). A PS task is refused by the engine rather
% than served FCFS, so it is refused here by name too. Element kinds, calls and
% precedences are declared for every model and are not marked.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = true;
reason = '';
if ~isa(model, 'LayeredNetwork')
    return
end
lsn = model.getStruct();
featSupported = SolverLDES.getLNFeatureSet();
featUsed = SolverFeatureSet;
nelem = lsn.nhosts + lsn.ntasks;
for idx = 1:nelem
    featUsed.setTrue(SchedStrategy.toFeature(lsn.sched(idx)));
    if idx > lsn.nhosts && lsn.sched(idx) == SchedStrategy.PS
        ok = false;
        reason = sprintf(['task ''%s'' is scheduled PS: a task holds threads and does not divide ' ...
            'them, so the LDES engine refuses it rather than serving it FCFS. Put PS on the ' ...
            'processor, or use FCFS, LCFS, SIRO, HOL or INF on the task.'], lsn.names{idx});
        return
    end
end
laws = {};
if isfield(lsn, 'hostdem') && iscell(lsn.hostdem)
    laws = [laws; lsn.hostdem(:)];
end
if isfield(lsn, 'think') && iscell(lsn.think)
    laws = [laws; lsn.think(:)];
end
for k = 1:numel(laws)
    law = laws{k};
    if isempty(law) || ~isobject(law) || ~ismethod(law, 'getFeatureName')
        continue
    end
    name = law.getFeatureName();
    if strcmp(name, 'Disabled')
        continue % an internal placeholder, as the Network recorder treats it
    end
    featUsed.setTrue(name);
end
[ok, reason] = SolverFeatureSet.supports(featSupported, featUsed);
if ~ok
    reason = sprintf(['The LDES layered engine does not serve this model. %s Use SolverLN, ' ...
        'or a discipline and a law from SolverLDES.getLNFeatureSet.'], strtrim(reason));
end
end
