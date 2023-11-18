function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn]=solver_ssa_analyzer_serial(sn, options, isHashed)
% [XN,UN,QN,RN,TN,CN]=SOLVER_SSA_ANALYZER_SERIAL(SN, OPTIONS)

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
schedid = sn.schedid;
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


options.samples = options.samples + 1;
[probSysState,StateSpaceAggr,arvRates,depRates,tranSysState,tranSync] = solver_ssa(sn, options);
%%
wset = 1:size(StateSpaceAggr,1);
for k=1:K
    refsf = sn.stationToStateful(sn.refstat(k));
    XN(k) = probSysState*depRates(wset,refsf,k);
end

for i=1:M
    isf = sn.stationToStateful(i);
    for k=1:K
        TN(i,k) = probSysState*depRates(wset,isf,k);
        QN(i,k) = probSysState*StateSpaceAggr(wset,(i-1)*K+k);
    end
    switch schedid(i)
        case SchedStrategy.ID_INF
            for k=1:K
                UN(i,k) = QN(i,k);
            end
        case {SchedStrategy.ID_PS, SchedStrategy.ID_DPS, SchedStrategy.ID_GPS}
            if isempty(sn.lldscaling) && isempty(sn.cdscaling)
                for k=1:K
                    if ~isempty(PH{i}{k})
                        UN(i,k) = probSysState*arvRates(wset,isf,k)*map_mean(PH{i}{k})/S(i);
                    end
                end
            else % lld/cd cases
                ind = sn.stationToNode(i);
                UN(i,1:K) = NaN;
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
            if isempty(sn.lldscaling) && isempty(sn.cdscaling)
                for k=1:K
                    if ~isempty(PH{i}{k})
                        UN(i,k) = probSysState*arvRates(wset,isf,k)*map_mean(PH{i}{k})/S(i);
                    end
                end
            else % lld/cd cases
                ind = sn.stationToNode(i);
                UN(i,1:K) = NaN;
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
    for i=1:M
        if TN(i,k)>0
            RN(i,k) = QN(i,k)./TN(i,k);
        else
            RN(i,k)=0;
        end
    end
    CN(k) = NK(k)./XN(k);
end
%%

% now update the routing probabilities in nodes with state-dependent routing
for k=1:K
    for isf=1:sn.nstateful
        if sn.nodetype(isf) == NodeType.Cache
            ind = sn.statefulToNode(isf);
            TNcache(isf,k) = probSysState*depRates(:,isf,k);
        end
    end
end

% updates cache actual hit and miss data
for k=1:K
    for isf=1:sn.nstateful
        if sn.nodetype(isf) == NodeType.Cache
            ind = sn.statefulToNode(isf);
            if length(sn.nodeparam{ind}.hitclass)>=k
                h = sn.nodeparam{ind}.hitclass(k);
                m = sn.nodeparam{ind}.missclass(k);
                if h==0 || m==0
                    sn.nodeparam{ind}.actualhitprob(k) = NaN;
                    sn.nodeparam{ind}.actualmissprob(k) = NaN;
                else
                    sn.nodeparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                    sn.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
                end
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
