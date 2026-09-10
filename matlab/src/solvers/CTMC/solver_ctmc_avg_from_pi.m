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
% clamp removes numerical residues, but an ME stationary vector is genuinely SIGNED, so clamping deletes real mass -- see solver_ctmc_analyzer
if ~(isfield(sn,'isph') && ~isempty(sn.isph) && ~all(sn.isph(:)))
    probSysState(probSysState<GlobalConstants.Zero) = 0;
end
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

% Column span of each STATION inside a StateSpace row; see the twin block in
% solver_ctmc_analyzer.m. sn.space is keyed by STATEFUL index and a stateful
% node need not be a station (a Cache is stateful and is not), so the running
% offset walks every stateful node and is read back through stationToStateful.
sfSpaceShift = zeros(1,sn.nstateful);
for isf=2:sn.nstateful
    sfSpaceShift(isf) = sfSpaceShift(isf-1) + size(sn.space{isf-1},2);
end
istSpaceShift = zeros(1,M);
istSpaceWidth = zeros(1,M);
for ist=1:M
    isf = sn.stationToStateful(ist);
    istSpaceShift(ist) = sfSpaceShift(isf);
    istSpaceWidth(ist) = size(sn.space{isf},2);
end

for k=1:K
    refsf = sn.stationToStateful(sn.refstat(k));
    XN(k) = probSysState*arvRates(wset,refsf,k);
end

% see _kb/06-solver-catalog.md (G-network signals) for rationale
inDropRegion = false(1,M);
if isfield(sn,'nregions') && sn.nregions > 0
    for f=1:sn.nregions
        if sn.regionrule(f) == DropStrategy.DROP
            memb = any(sn.region{f} ~= -1, 2); % stations constrained by region f
            memb = memb(:)';
            inDropRegion(1:min(M,numel(memb))) = inDropRegion(1:min(M,numel(memb))) | memb(1:min(M,numel(memb)));
        end
    end
end

for ist=1:M
    isf = sn.stationToStateful(ist);
    ind = sn.stationToNode(ist);
    isSource = sn.nodetype(ind) == NodeType.Source;
    for k=1:K
        TN(ist,k) = probSysState*depRates(wset,isf,k);
        if isSource
            % State.toMarginal encodes an EXT station as nir = Inf, an infinite
            % reservoir, which is a statement about the state space and not a
            % queue length. Reading it as one gave Q = Inf at the Source.
            QN(ist,k) = 0;
        else
            QN(ist,k) = probSysState*StateSpaceAggr(wset,(ist-1)*K+k);
        end
    end
    if ~isSource
        % see _kb/06-solver-catalog.md (G-network signals) for rationale
        % A class that can be DROPPED here must be measured on the CARRIED rate
        % alone: the offered rate counts arrivals that never entered service, so
        % max(UNarv,UNdep) would report the offered load as utilization (an
        % M/M/1/4 with lambda=0.6 gave 0.6 against the true 1-p0=0.566). This
        % guard is present in solver_ctmc_analyzer and was missing here, so the
        % two implementations disagreed on every lossy station -- the same drift
        % recorded for this file's lld/cd branch on 2026-07-17.
        canDropClass = isinf(sn.njobs(:)') & (isfinite(sn.cap(ist)) | isfinite(sn.classcap(ist,:)) | inDropRegion(ist));
        signalLoss = ctmc_signal_lossy(sn, arvRates, probSysState, wset, isf);
        canDropClass = canDropClass | signalLoss;
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
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            if canDropClass(k)
                                UN(ist,k) = UNdep_ik;
                            else
                                UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                                UN(ist,k) = max(UNarv_ik,UNdep_ik);
                            end
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
                        [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
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
                    [~,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
                    for k=1:K
                        UN(ist,k) = UN(ist,k) + probSysState(st)*sir(k)/S(ist);
                    end
                end
            otherwise
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            if canDropClass(k)
                                UN(ist,k) = UNdep_ik;
                            else
                                UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                                UN(ist,k) = max(UNarv_ik,UNdep_ik);
                            end
                        end
                    end
                else % lld/cd/ljd cases
                    ind = sn.stationToNode(ist);
                    % the load-dependent station's capacity is its PEAK scaling,
                    % not its server count, so normalize by ceff not S
                    ceff = S(ist);
                    if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                        ceff = max(ceff, max(sn.lldscaling(ist,:)));
                    end
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
                        if ni>0
                            lldnow = 1;
                            if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                                lldnow = sn.lldscaling(ist, min(max(sum(ni),1), size(sn.lldscaling,2)));
                            end
                            sirtot = sum(sir);
                            for k=1:K
                                if sirtot > 0
                                    UN(ist,k) = UN(ist,k) + probSysState(st)*(sir(k)/sirtot)*lldnow/ceff;
                                end
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

% Class-, joint- and global-dependence utilization normalization Util=T*S/peak,
% using the declared sn.cdscalingpeak, sn.jdscalingpeak and sn.gdscalingpeak;
% see _kb/06-solver-catalog.md (Utilization conventions) for rationale.
hasgd = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
if (~isempty(sn.cdscaling) || ~isempty(sn.jdscaling) || hasgd) && ~any(isinf(sn.njobs))
    for ist=1:M
        hascd = ~isempty(sn.cdscaling) && length(sn.cdscaling) >= ist && ~isempty(sn.cdscaling{ist});
        hasjd = ~isempty(sn.jdscaling) && length(sn.jdscaling) >= ist && ~isempty(sn.jdscaling{ist});
        if ~hascd && ~hasjd && ~hasgd
            continue
        end
        for k=1:K
            % beta_r(n), eta_i(n) and phi(n) scale the SAME rate, so the peaks multiply
            bmax = 1;
            if hascd
                bmax = bmax * sn.cdscalingpeak(ist,k);
            end
            if hasjd
                bmax = bmax * sn.jdscalingpeak(ist,k);
            end
            if hasgd
                bmax = bmax * sn.gdscalingpeak(ist,k);
            end
            if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
            else
                UN(ist,k) = 0;
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
