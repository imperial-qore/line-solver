function isReneging = mam_has_reneging_patience(sn)
% ISRENEGING = MAM_HAS_RENEGING_PATIENCE(SN)
%
% True when some station-class pair declares RENEGING impatience together with
% a patience distribution, i.e. when solver_mam_retrial has a MAP/M/s+G model
% to solve rather than a BMAP/PH/N/N retrial one.
%
% Lifted out of solver_mam_analyzer.m so that SolverMAM.supportsModelMethod can
% ask the same question: the gate decides whether to OFFER the 'retrial' method
% and the analyzer decides whether to run it, and two copies of the test are
% how the two come to disagree.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

isReneging = false;

% ImpatienceType (RENEGING, BALKING) per station-class; absent on a model that
% declares no impatience at all.
if ~isfield(sn, 'impatienceClass') || isempty(sn.impatienceClass)
    return;
end

% The patience law itself. impatienceClass is set from the station's impatience
% type independently of whether a distribution was configured, so both fields
% have to be present before the pair counts as a reneging model.
if ~isfield(sn, 'patienceProc') || isempty(sn.patienceProc)
    return;
end

for ist = 1:sn.nstations
    for r = 1:sn.nclasses
        if sn.impatienceClass(ist, r) == ImpatienceType.RENEGING
            if ~isempty(sn.patienceProc{ist, r})
                isReneging = true;
                return;
            end
        end
    end
end
end
