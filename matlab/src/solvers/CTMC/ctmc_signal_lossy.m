function lossy = ctmc_signal_lossy(sn, arvRates, probSysState, wset, isf)
% LOSSY = CTMC_SIGNAL_LOSSY(SN, ARVRATES, PROBSYSSTATE, WSET, ISF)
%
% Classes whose jobs can be annihilated by a G-network signal at the stateful
% node ISF, as a 1 x nclasses logical row. Such a job leaves the station
% without a service completion, so the arrival-based utilization estimator
% (offered load) is invalid there and only the departure-based (carried load)
% estimator UN = T * E[S] / c is meaningful.
%
% A signal class R is active at ISF when its stationary arrival rate there is
% positive. A targeted signal (SN.SIGNALTARGET(R) >= 1) only removes that
% class; an untargeted one is class-agnostic and removes any non-signal class,
% matching State.afterEventStationSignal, MAM and LDES.
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

K = sn.nclasses;
lossy = false(1,K);
if ~isfield(sn,'issignal') || isempty(sn.issignal) || ~any(sn.issignal)
    return
end
if isempty(arvRates) || isempty(wset)
    return
end
issignal = logical(sn.issignal(:)');
for r=1:K
    if ~issignal(r)
        continue
    end
    if probSysState*arvRates(wset,isf,r) <= 0
        continue
    end
    tgt = -1;
    if isfield(sn,'signaltarget') && numel(sn.signaltarget) >= r
        tgt = sn.signaltarget(r);
    end
    if tgt >= 1 && tgt <= K
        lossy(tgt) = true;
    else
        lossy(~issignal) = true;
    end
end
end
