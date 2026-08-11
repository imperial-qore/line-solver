function [QN,UN,RN,TN,CN,XN] = solver_ctmc_avg_from_pi(sn, pivec, StateSpace, StateSpaceAggr, arvRates, depRates, options)
% SOLVER_CTMC_AVG_FROM_PI  Map a state distribution to mean performance metrics.
%
% [QN,UN,RN,TN,CN,XN] = SOLVER_CTMC_AVG_FROM_PI(SN, PIVEC, STATESPACE,
%   STATESPACEAGGR, ARVRATES, DEPRATES)
%
% Given an arbitrary probability vector PIVEC over the enumerated CTMC state
% space of SN (rows of STATESPACE / STATESPACEAGGR), returns the per-(station,
% class) mean queue length QN, utilization UN, response time RN, throughput TN,
% system response time CN and system throughput XN. The discipline-aware mapping
% is identical to the steady-state reduction performed by solver_ctmc_analyzer;
% it is factored here so that callers holding their own distribution (e.g. the
% SolverENV state-vector analyzer, which time-averages a transient distribution)
% can reuse it without re-solving for the stationary vector.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
S = sn.nservers;
NK = sn.njobs';
sched = sn.sched;
PH = sn.proc;

probSysState = pivec(:)';
probSysState(probSysState<GlobalConstants.Zero) = 0;
if sum(probSysState) > 0
    probSysState = probSysState/sum(probSysState);
end
wset = 1:size(StateSpace,1);

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);

istSpaceShift = zeros(1,M);
for ist=1:M
    if ist==1
        istSpaceShift(ist) = 0;
    else
        istSpaceShift(ist) = istSpaceShift(ist-1) + size(sn.space{ist-1},2);
    end
end

for k=1:K
    refsf = sn.stationToStateful(sn.refstat(k));
    XN(k) = probSysState*arvRates(wset,refsf,k);
end

for ist=1:M
    isf = sn.stationToStateful(ist);
    ind = sn.stationToNode(ist);
    for k=1:K
        TN(ist,k) = probSysState*depRates(wset,isf,k);
        QN(ist,k) = probSysState*StateSpaceAggr(wset,(ist-1)*K+k);
    end
    if sn.nodetype(ind) ~= NodeType.Source
        % see _kb/06-solver-catalog.md (G-network signals) for rationale
        signalLoss = ctmc_signal_lossy(sn, arvRates, probSysState, wset, isf);
        switch sched(ist)
            case SchedStrategy.INF
                for k=1:K
                    UN(ist,k) = QN(ist,k);
                end
            case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS}
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                            UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            UN(ist,k) = signalLoss(k)*UNdep_ik + (1-signalLoss(k))*max(UNarv_ik,UNdep_ik);
                        end
                    end
                else % lld/cd/ljd cases
                    % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                    ind = sn.stationToNode(ist);
                    ceff = S(ist);
                    if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                        ceff = max(ceff, max(sn.lldscaling(ist,:)));
                    end
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
                        if ni>0
                            lldnow = 1;
                            if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                                lldnow = sn.lldscaling(ist, min(max(sum(ni),1), size(sn.lldscaling,2)));
                            end
                            for k=1:K
                                UN(ist,k) = UN(ist,k) + probSysState(st)*nir(k)*sn.schedparam(ist,k)/(nir*sn.schedparam(ist,:)')*lldnow/ceff;
                            end
                        end
                    end
                end
            case SchedStrategy.PAS
                % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                ind = sn.stationToNode(ist);
                UN(ist,1:K) = 0;
                for st = wset
                    [~,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
                    for k=1:K
                        UN(ist,k) = UN(ist,k) + probSysState(st)*sir(k)/S(ist);
                    end
                end
            otherwise
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                            UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            UN(ist,k) = signalLoss(k)*UNdep_ik + (1-signalLoss(k))*max(UNarv_ik,UNdep_ik);
                        end
                    end
                else % lld/cd/ljd cases
                    ind = sn.stationToNode(ist);
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
                        if ni>0
                            for k=1:K
                                UN(ist,k) = UN(ist,k) + probSysState(st)*sir(k)/S(ist);
                            end
                        end
                    end
                end
        end
        % see _kb/06-solver-catalog.md (G-network signals) for rationale
        if any(signalLoss) && sched(ist) ~= SchedStrategy.INF && ...
                isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
            UNb = ctmc_signal_busy(sn, ind, ist, sched(ist), S(ist), StateSpace, istSpaceShift, wset, probSysState);
            UN(ist,signalLoss) = UNb(signalLoss);
        end
    end
end

% see _kb/06-solver-catalog.md (Utilization conventions) for rationale
if ~isempty(sn.cdscaling) && ~any(isinf(sn.njobs))
    for ist=1:M
        if length(sn.cdscaling) >= ist && ~isempty(sn.cdscaling{ist})
            for k=1:K
                bmax = sn.cdscalingpeak(ist,k);
                if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                    UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
                else
                    UN(ist,k) = 0;
                end
            end
        end
    end
end

% Joint-dependence (non-product-form) utilization normalization: Util=T*S/peak
% using the declared sn.jdscalingpeak, mirroring the class-dependence block.
if ~isempty(sn.jdscaling) && ~any(isinf(sn.njobs))
    for ist=1:M
        if length(sn.jdscaling) >= ist && ~isempty(sn.jdscaling{ist})
            for k=1:K
                bmax = sn.jdscalingpeak(ist,k);
                if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                    UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
                else
                    UN(ist,k) = 0;
                end
            end
        end
    end
end


for k=1:K
    for ist=1:M
        if TN(ist,k)>0
            RN(ist,k) = QN(ist,k)./TN(ist,k);
        else
            RN(ist,k)=0;
        end
    end
    CN(k) = NK(k)./XN(k);
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;
end
