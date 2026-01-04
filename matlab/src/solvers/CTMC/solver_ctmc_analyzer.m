function [QN,UN,RN,TN,CN,XN,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname,sncopy] = solver_ctmc_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME,sn] = SOLVER_CTMC_ANALYZER(sn, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%if options.remote
%    sn.rtfun = {};
%    sn.lst = {};
%    qn_json = jsonencode(sn);
%    sn = NetworkStruct.fromJSON(qn_json)
%return
%end

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
S = sn.nservers;
NK = sn.njobs';  % initial population per class
sched = sn.sched;

Tstart = tic;
PH = sn.proc;

line_debug('CTMC analyzer starting: nstations=%d, nclasses=%d, njobs=%s', M, K, mat2str(NK));

% Note: hide_immediate now selectively preserves Cache immediate transitions
% in solver_ctmc.m, so we no longer need to disable it entirely for Cache nodes

line_debug('Building state space and infinitesimal generator via solver_ctmc');
[InfGen,StateSpace,StateSpaceAggr,EventFiltration,arvRates,depRates,sn] = solver_ctmc(sn, options); % sn is updated with the state space

% if the initial state does not reflect the final size of the state
% vectors, attempt to correct it
for isf=1:sn.nstateful
    if size(sn.state{isf},2) < size(sn.space{isf},2)
        sn.state{isf} = [zeros(1,size(sn.space{isf},2)-size(sn.state{isf},2)),sn.state{isf}];
    end
end
sncopy = sn;

if options.keep
    line_debug('Saving CTMC data to file (options.keep=true)');
    fname = lineTempName;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    line_printf('CTMC infinitesimal generator and state space saved in: ');
    line_printf(strrep(sprintf('%s.mat\n',fname),'\','\\'))
else
    fname = '';
end

wset = 1:length(InfGen);

line_debug('State space built: %d states, solving CTMC', length(InfGen));

use_ctmc_solve_stable = true;
if use_ctmc_solve_stable
    % stable version
    % Note: solver_ctmc now selectively preserves Cache immediate transitions
    % to enable hit/miss rate computation while hiding other immediate transitions
    [probSysState, ~, nConnComp, connComp] = ctmc_solve(InfGen, options);

    if nConnComp > 1
        line_debug('CTMC is reducible: %d connected components', nConnComp);
        % the matrix was reducible
        initState = matchrow(StateSpace, cell2mat(sn.state'));
        % determine the weakly connected component associated to the initial state
        wset = find(connComp == connComp(initState));
        line_debug('Using component %d with %d states (from initial state)', connComp(initState), length(wset));
        probSysState = ctmc_solve(InfGen(wset, wset), options);
        InfGen = InfGen(wset, wset);
        StateSpace = StateSpace(wset,:);
    else
        line_debug('CTMC is irreducible, using full state space');
    end
else
    % development version

    % we now find the initial state and then solver the CTMC allowing for the
    % case where it is reducible
    initState = matchrow(StateSpace, cell2mat(sn.state'));
    pi0 = zeros(1,length(InfGen)); pi0(initState) = 1.0;
    [pi,pis,~,scc,~] = ctmc_solve_reducible(InfGen, pi0, options);

    if size(pis,1)==1
        probSysState = pi;
    else
        wset = scc == scc(initState);
        InfGen = InfGen(wset, wset);
        StateSpace = StateSpace(wset,:);
        probSysState = pis(scc(initState),scc == scc(initState));
    end
end
probSysState(probSysState<GlobalConstants.Zero)=0;
probSysState = probSysState/sum(probSysState);

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
        switch sched(ist)
            case SchedStrategy.INF
                for k=1:K
                    UN(ist,k) = QN(ist,k);
                end
            case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS}
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && ~sn_has_joint_dependence(sn)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % There are cases where due to remove of
                            % immediate transitions or due to cutoff the
                            % utilization estimator based on arrivals or
                            % departures can be under-estimated. E.g.:
                            % UNarv_ik doesn't work well with example_stateDependentRouting_3
                            % UNdep_ik doesn't work well with test_OQN_JMT_6
                            % Therefore, we take the maximum of the two
                            % Note: the two estimators are normally
                            % identical.
                            UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            UN(ist,k) = max(UNarv_ik,UNdep_ik);
                        end
                    end
                else % lld/cd/ljd cases
                    ind = sn.stationToNode(ist);
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
                        if ni>0
                            for k=1:K
                                UN(ist,k) = UN(ist,k) + probSysState(st)*nir(k)*sn.schedparam(ist,k)/(nir*sn.schedparam(ist,:)');
                            end
                        end
                    end
                end
            otherwise
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && ~sn_has_joint_dependence(sn)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % There are cases where due to remove of
                            % immediate transitions or due to cutoff the
                            % utilization estimator based on arrivals or
                            % departures can be under-estimated. E.g.:
                            % UNarv_ik doesn't work well with example_stateDependentRouting_3
                            % UNdep_ik doesn't work well with test_OQN_JMT_6
                            % Therefore, we take the maximum of the two
                            % Note: the two estimators are normally
                            % identical.
                            UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            UN(ist,k) = max(UNarv_ik,UNdep_ik);
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

runtime = toc(Tstart);

% now update the routing probabilities in nodes with state-dependent routing
TNcache = zeros(sn.nstateful,K);
for k=1:K
    for isf=1:sn.nstateful
        ind = sncopy.statefulToNode(isf);
        if sncopy.nodetype(ind) == NodeType.Cache
            TNcache(isf,k) = probSysState*depRates(wset,isf,k);
        end
    end
end

% updates cache actual hit and miss data
for k=1:K
    for isf=1:sncopy.nstateful
        ind = sncopy.statefulToNode(isf);
        if sncopy.nodetype(ind) == NodeType.Cache
            if length(sncopy.nodeparam{ind}.hitclass)>=k
                h = sncopy.nodeparam{ind}.hitclass(k);
                m = sncopy.nodeparam{ind}.missclass(k);
                if h> 0 && m > 0
                    sncopy.nodeparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                    sncopy.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
                end
            end
        end
    end
end
end
