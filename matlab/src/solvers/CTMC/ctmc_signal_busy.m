function UNb = ctmc_signal_busy(sn, ind, ist, schedIst, SIst, StateSpace, istSpaceShift, wset, probSysState)
% UNB = CTMC_SIGNAL_BUSY(SN, IND, IST, SCHEDIST, SIST, STATESPACE,
%   ISTSPACESHIFT, WSET, PROBSYSSTATE)
%
% Exact per-class busy-server fraction at station IST, read off the enumerated
% state space as a 1 x nclasses row.
%
% For a class that a G-network signal can annihilate, the departure-based
% estimator T*E[S]/c is exact only under exponential service: a job destroyed
% mid-service leaves behind busy time with no completion, so with phase-type
% service T*E[S]/c under-counts (M/Er2/1 with lambda+=0.5, lambda-=0.4 gives
% 0.34941 against a true 0.37696). The in-service occupancy below is exact for
% any service process.
%
% PS-like disciplines share the servers among all resident jobs, so class k
% gets the weighted share n_k w_k / sum_j n_j w_j of the busy servers; the
% remaining disciplines expose the in-service indicator directly through
% State.toMarginal.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
UNb = zeros(1,K);
isPS = any(schedIst == [SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS]);
cols = (istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2));
for st = wset
    if probSysState(st) == 0
        continue
    end
    [ni,nir,sir] = State.toMarginal(sn, ind, StateSpace(st,cols));
    if ni <= 0
        continue
    end
    if isPS
        w = sn.schedparam(ist,:);
        wtot = nir(:)'*w(:);
        if wtot > 0
            for k=1:K
                UNb(k) = UNb(k) + probSysState(st)*(nir(k)*w(k)/wtot)*min(sum(ni),SIst)/SIst;
            end
        end
    else
        for k=1:K
            UNb(k) = UNb(k) + probSysState(st)*sir(k)/SIst;
        end
    end
end
end
