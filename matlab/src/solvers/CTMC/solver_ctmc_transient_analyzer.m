function [t,pit,QNt,UNt,RNt,TNt,CNt,XNt,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname] = solver_ctmc_transient_analyzer(sn, options)
% [T,PIT,QNT,UNT,RNT,TNT,CNT,XNT,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME] = SOLVER_CTMC_TRANSIENT_ANALYZER(QN, OPTIONS)
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

RNt=[]; CNt=[];  XNt=[];

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
fname = '';
Tstart = tic;
S = sn.nservers;
schedid = sn.schedid;
PH = sn.proc;

[InfGen,StateSpace,StateSpaceAggr,EventFiltration,~,depRates,sn] = solver_ctmc(sn, options); % sn is updated with the state space

if options.keep
    fname = lineTempName;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    line_printf('\nCTMC infinitesimal generator and state space saved in: ');
    line_printf([fname, '.mat'])
end

state = [];
for i=1:sn.nnodes
    if sn.isstateful(i)
        isf = sn.nodeToStateful(i);
        state = [state,zeros(1,size(sn.space{isf},2)-length(sn.state{isf})),sn.state{isf}];
    end
end
pi0 = zeros(1,length(InfGen));

state0 = matchrow(StateSpace, state);
if state0 == -1
    line_error(mfilename,'Initial state not contained in the state space.');   
%     state0 = matchrow(StateSpace, round(state));
%     state = round(state);
%     if state0 == -1
%         line_error(mfilename,'Cannot recover - CTMC stopping');
%     end
end
pi0(state0) = 1; % find initial state and set it to probability 1

%if options.timespan(1) == options.timespan(2)
%    pit = ctmc_uniformization(pi0,Q,options.timespan(1));
%    t = options.timespan(1);
%else
[pit,t] = ctmc_transient(InfGen,pi0,options.timespan(1),options.timespan(2),options.stiff);
%end
pit(pit<GlobalConstants.Zero)=0;

QNt = cell(M,K);
UNt = cell(M,K);
%XNt = cell(1,K);
TNt = cell(M,K);

if t(1) == 0
    t(1) = GlobalConstants.Zero;
end
for k=1:K
    %    XNt(k) = pi*arvRates(:,sn.refstat(k),k);
    for i=1:M        
        %occupancy_t = cumsum(pit.*[0;diff(t)],1)./t;        
        occupancy_t = pit;
        TNt{i,k} = occupancy_t*depRates(:,i,k);
        qlenAt_t = pit*StateSpaceAggr(:,(i-1)*K+k);
        %QNt{i,k} = cumsum(qlenAt_t.*[0;diff(t)])./t;
        QNt{i,k} = qlenAt_t;
        switch schedid(i)
            case SchedStrategy.ID_INF
                UNt{i,k} = QNt{i,k};
            case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT, SchedStrategy.ID_SJF}
                if ~isempty(PH{i}{k})
                    UNt{i,k} = occupancy_t*min(StateSpaceAggr(:,(i-1)*K+k),S(i))/S(i);
                end
            case SchedStrategy.ID_PS
                uik = min(StateSpaceAggr(:,(i-1)*K+k),S(i)) .* StateSpaceAggr(:,(i-1)*K+k) ./ sum(StateSpaceAggr(:,((i-1)*K+1):(i*K)),2);
                uik(isnan(uik))=0;
                utilAt_t = pit * uik / S(i);
                %UNt{i,k} = cumsum(utilAt_t.*[0;diff(t)])./t;
                UNt{i,k} = utilAt_t;
            case SchedStrategy.ID_DPS
                w = sn.schedparam(i,:);
                nik = S(i) * w(k) * StateSpaceAggr(:,(i-1)*K+k) ./ sum(repmat(w,size(StateSpaceAggr,1),1).*StateSpaceAggr(:,((i-1)*K+1):(i*K)),2);
                nik(isnan(nik))=0;
                UNt{i,k} = occupancy_t*nik;
            otherwise
                if ~isempty(PH{i}{k})
                    ind = sn.stationToNode(i);
                    line_warning(mfilename,'Transient utilization not support yet for station %s, returning an approximation.\n',sn.nodenames{ind});
                    UNt{i,k} = occupancy_t*min(StateSpaceAggr(:,(i-1)*K+k),S(i))/S(i);
                end
        end
    end
end
runtime = toc(Tstart);

%if options.verbose
%    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
