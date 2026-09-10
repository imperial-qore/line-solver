function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN,tranStartTag,tranPreemptTag]=solver_ssa_analyzer_serial(sn, init_state, options, isHashed)
% [XN,UN,QN,RN,TN,CN,TRANSYSSTATE,TRANSYNC,SN,STARTN,PREEMPTN,TRANSTARTTAG,TRANPREEMPTTAG]=SOLVER_SSA_ANALYZER_SERIAL(SN, OPTIONS)
%
% STARTN and PREEMPTN are the (station x class) service-start and preemption
% rates, estimated exactly as TN is: the time average of the rate the enabled
% transitions carry. TRANSTARTTAG/TRANPREEMPTTAG are the per-step tags of the
% sampled path, used by @SolverSSA/sample.m to emit START and PREEMPT events.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes

% istSpaceShift = zeros(1,M);
% for i=1:M
%     if i==1
%         istSpaceShift(i) = 0;
%     else
%         istSpaceShift(i) = istSpaceShift(i-1) + size(sn.space{i-1},2);
%     end
% end

S = sn.nservers;
NK = sn.njobs';  % initial population per class
sched = sn.sched;
Tstart = tic;

PH = sn.proc;
tranSysState = [];
tranSync = [];

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);


if isfield(options,'config') && isfield(options.config,'eventcache')
    eventCache = EventCache.create(options.config.eventcache, sn);
else
    eventCache = EventCache.create(false, sn);
end

% see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
userCap = sn.cap;
userClasscap = sn.classcap;
[probSysState,StateSpaceAggr,arvRates,depRates,tranSysState,tranSync,~,dlyRates,startRates,preemptRates,tranStartTag,tranPreemptTag] = solver_ssa(sn, init_state, options, eventCache);
%%
wset = 1:size(StateSpaceAggr,1);
StartN = zeros(M,K);
PreemptN = zeros(M,K);
for k=1:K
    refsf = sn.stationToStateful(sn.refstat(k));
    XN(k) = probSysState*depRates(wset,refsf,k);
end

% Global (Whittle) rate scaling: like cd/jd it rescales the service rate, so
% the peak it declares is what normalizes the utilization.
hasGd = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
for ist=1:M
    isf = sn.stationToStateful(ist);
    for k=1:K
        TN(ist,k) = probSysState*depRates(wset,isf,k);
        QN(ist,k) = probSysState*StateSpaceAggr(wset,(ist-1)*K+k);
        % same time average as TN, over the derived tag rates
        StartN(ist,k) = probSysState*startRates(wset,isf,k);
        PreemptN(ist,k) = probSysState*preemptRates(wset,isf,k);
    end
    % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
    canDropClass = isinf(sn.njobs(:)') & (isfinite(userCap(ist)) | isfinite(userClasscap(ist,:)));
    switch sched(ist)
        case SchedStrategy.INF
            for k=1:K
                UN(ist,k) = QN(ist,k);
            end
        case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
              SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
            if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling) && ~hasGd
                for k=1:K
                    if ~isempty(PH{ist}{k})
                        if canDropClass(k)
                            UN(ist,k) = TN(ist,k)/sn.rates(ist,k)/S(ist);
                        else
                            UN(ist,k) = probSysState*arvRates(wset,isf,k)/sn.rates(ist,k)/S(ist);
                        end
                    end
                end
            else % lld/cd/ljd cases
                ind = sn.stationToNode(ist);
                % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
                isCd = ~isempty(sn.cdscaling) && ist <= numel(sn.cdscaling) && ~isempty(sn.cdscaling{ist});
                isJd = ~isempty(sn.jdscaling) && ist <= numel(sn.jdscaling) && ~isempty(sn.jdscaling{ist});
                isGd = hasGd;
                ceff = S(ist);
                if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                    ceff = max(ceff, max(sn.lldscaling(ist,:)));
                end
                for k=1:K
                    if ~isempty(PH{ist}{k})
                        if isCd || isJd || isGd
                            % effective peak = product of declared cd and jd peaks
                            cdiv = 1;
                            if isCd, cdiv = cdiv * sn.cdscalingpeak(ist,k); end
                            if isJd, cdiv = cdiv * sn.jdscalingpeak(ist,k); end
                            if isGd, cdiv = cdiv * sn.gdscalingpeak(ist,k); end
                        else
                            cdiv = ceff;
                        end
                        if cdiv > 0
                            UN(ist,k) = TN(ist,k)*map_mean(PH{ist}{k})/cdiv;
                        else
                            UN(ist,k) = 0;
                        end
                    end
                end
                %               UN(i,1:K) = 0;
                %                 for st = wset
                %                     [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(i)+1):(istSpaceShift(i)+size(sn.space{i},2))));
                %                     if ni>0
                %                         for k=1:K
                %                             UN(i,k) = UN(i,k) + probSysState(st)*nir(k)*sn.schedparam(i,k)/(nir*sn.schedparam(i,:)');
                %                         end
                %                     end
                %                 end
            end
        otherwise
            if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling) && ~hasGd
                for k=1:K
                    if ~isempty(PH{ist}{k})
                        if canDropClass(k)
                            UN(ist,k) = TN(ist,k)*map_mean(PH{ist}{k})/S(ist);
                        else
                            UN(ist,k) = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                        end
                    end
                end
            else % lld/cd/ljd cases
                ind = sn.stationToNode(ist);
                % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
                isCd = ~isempty(sn.cdscaling) && ist <= numel(sn.cdscaling) && ~isempty(sn.cdscaling{ist});
                isJd = ~isempty(sn.jdscaling) && ist <= numel(sn.jdscaling) && ~isempty(sn.jdscaling{ist});
                isGd = hasGd;
                ceff = S(ist);
                if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                    ceff = max(ceff, max(sn.lldscaling(ist,:)));
                end
                for k=1:K
                    if ~isempty(PH{ist}{k})
                        if isCd || isJd || isGd
                            % effective peak = product of declared cd and jd peaks
                            cdiv = 1;
                            if isCd, cdiv = cdiv * sn.cdscalingpeak(ist,k); end
                            if isJd, cdiv = cdiv * sn.jdscalingpeak(ist,k); end
                            if isGd, cdiv = cdiv * sn.gdscalingpeak(ist,k); end
                        else
                            cdiv = ceff;
                        end
                        if cdiv > 0
                            UN(ist,k) = TN(ist,k)*map_mean(PH{ist}{k})/cdiv;
                        else
                            UN(ist,k) = 0;
                        end
                    end
                end
                %               UN(i,1:K) = 0;
                %                 for st = wset
                %                     [ni,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(i)+1):(istSpaceShift(i)+size(sn.space{i},2))));
                %                     if ni>0
                %                         for k=1:K
                %                             UN(i,k) = UN(i,k) + probSysState(st)*sir(k)/S(i);
                %                         end
                %                     end
                %                 end
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
%%

% now update the routing probabilities in nodes with state-dependent routing
TNcache = zeros(sn.nstateful, K);
XNcache = zeros(sn.nstateful, K);
DNcache = zeros(sn.nstateful, K);
for k=1:K
    for isf=1:sn.nstateful
        if sn.nodetype(isf) == NodeType.Cache
            TNcache(isf,k) = probSysState*depRates(:,isf,k);
            XNcache(isf,k) = probSysState*arvRates(:,isf,k);
            if ~isempty(dlyRates)
                DNcache(isf,k) = probSysState*dlyRates(:,isf,k);
            end
        end
    end
end

% see _kb/09-ldes-and-cache.md for rationale (SSA cache hit/miss accounting)
retrievalLatencyWarned = false;
for k=1:K
    for isf=1:sn.nstateful
        if sn.nodetype(isf) == NodeType.Cache
            ind = sn.statefulToNode(isf);
            np = sn.nodeparam{ind};
            if length(np.hitclass)>=k
                h = np.hitclass(k);
                m = np.missclass(k);
                if h>0 && m>0
                    tot = sum(TNcache(isf,[h,m]));
                    % The hit-class departure rate already contains the released
                    % delayed hits: a request merged onto an in-flight fetch is
                    % released in the hit class when the fetch completes. Carve
                    % them out rather than adding a fourth share.
                    dly = DNcache(isf,k);
                    sn.nodeparam{ind}.actualhitprob(k) = max(TNcache(isf,h)-dly,0)/tot;
                    sn.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/tot;
                    if dly > 0
                        sn.nodeparam{ind}.actualdelayedhitprob(k) = dly/tot;
                    end
                else
                    sn.nodeparam{ind}.actualhitprob(k) = NaN;
                    sn.nodeparam{ind}.actualmissprob(k) = NaN;
                end

                % see _kb/09-ldes-and-cache.md for rationale (SSA retrieval latency NaN)
                expectedLatency = NaN;
                if isfield(np, 'retrievalSystemQueueIndices') ...
                        && isKey(np.retrievalSystemQueueIndices, int32(k-1)) ...
                        && ~isempty(np.retrievalSystemQueueIndices{int32(k-1)})
                    if ~retrievalLatencyWarned
                        line_warning(mfilename, 'Retrieval-system expected latency is not currently implemented; reporting NaN.');
                        retrievalLatencyWarned = true;
                    end
                end
                sn.nodeparam{ind}.actualresidt(k) = expectedLatency;
            end
        end
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;
end
