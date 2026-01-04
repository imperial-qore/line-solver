function [t,pit,QNt,UNt,RNt,TNt,CNt,XNt,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname] = solver_ctmc_transient_analyzer(sn, options)
% [T,PIT,QNT,UNT,RNT,TNT,CNT,XNT,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME] = SOLVER_CTMC_TRANSIENT_ANALYZER(QN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

RNt=[]; CNt=[];  XNt=[];

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
fname = '';
Tstart = tic;
S = sn.nservers;
sched = sn.sched;
PH = sn.proc;

line_debug('CTMC transient analyzer starting: nstations=%d, nclasses=%d', M, K);

[InfGen,StateSpace,StateSpaceAggr,EventFiltration,~,depRates,sn] = solver_ctmc(sn, options); % sn is updated with the state space

if options.keep
    fname = lineTempName;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    line_printf('\nCTMC infinitesimal generator and state space saved in: ');
    line_printf([fname, '.mat'])
end

state = [];
for ist=1:sn.nnodes
    if sn.isstateful(ist)
        isf = sn.nodeToStateful(ist);
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
[pit,t] = ctmc_transient(InfGen,pi0,options.timespan(1),options.timespan(2),options.stiff,[],options.timestep);
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
    for ist=1:M
        %occupancy_t = cumsum(pit.*[0;diff(t)],1)./t;
        occupancy_t = pit;
        TNt{ist,k} = occupancy_t*depRates(:,ist,k);
        qlenAt_t = pit*StateSpaceAggr(:,(ist-1)*K+k);
        %QNt{i,k} = cumsum(qlenAt_t.*[0;diff(t)])./t;
        QNt{ist,k} = qlenAt_t;
        switch sched(ist)
            case SchedStrategy.INF
                UNt{ist,k} = QNt{ist,k};
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.SJF}
                if ~isempty(PH{ist}{k})
                    UNt{ist,k} = occupancy_t*min(StateSpaceAggr(:,(ist-1)*K+k),S(ist))/S(ist);
                end
            case SchedStrategy.PS
                uik = min(StateSpaceAggr(:,(ist-1)*K+k),S(ist)) .* StateSpaceAggr(:,(ist-1)*K+k) ./ sum(StateSpaceAggr(:,((ist-1)*K+1):(ist*K)),2);
                uik(isnan(uik))=0;
                utilAt_t = pit * uik / S(ist);
                %UNt{i,k} = cumsum(utilAt_t.*[0;diff(t)])./t;
                UNt{ist,k} = utilAt_t;
            case SchedStrategy.DPS
                w = sn.schedparam(ist,:);
                nik = S(ist) * w(k) * StateSpaceAggr(:,(ist-1)*K+k) ./ sum(repmat(w,size(StateSpaceAggr,1),1).*StateSpaceAggr(:,((ist-1)*K+1):(ist*K)),2);
                nik(isnan(nik))=0;
                UNt{ist,k} = occupancy_t*nik;
            otherwise
                if ~isempty(PH{ist}{k})
                    ind = sn.stationToNode(ist);
                    line_warning(mfilename,'Transient utilization not support yet for station %s, returning an approximation.\n',sn.nodenames{ind});
                    UNt{ist,k} = occupancy_t*min(StateSpaceAggr(:,(ist-1)*K+k),S(ist))/S(ist);
                end
        end
    end
end
runtime = toc(Tstart);

%if options.verbose
%    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
