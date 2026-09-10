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
    state0 = matchrow(StateSpace, round(state));
    state = round(state);
    if state0 == -1
        line_error(mfilename,'Initial state not contained in the state space.');
    end
end
pi0(state0) = 1; % find initial state and set it to probability 1

% see _kb/06-solver-catalog.md (CTMC section, time-inhomogeneous generator) for rationale
rate_sched = [];
if isfield(options,'config') && isfield(options.config,'rate_sched') && ~isempty(options.config.rate_sched)
    rate_sched = options.config.rate_sched;
end
% see _kb/06-solver-catalog.md (CTMC section, transient methods) for rationale
transient_method = 'ode';
if isfield(options,'config') && isfield(options.config,'transient_method') && ~isempty(options.config.transient_method)
    transient_method = lower(options.config.transient_method);
end
if isempty(rate_sched)
    switch transient_method
        case 'fau'
            [pit,t] = local_ctmc_fau(InfGen,pi0,options);
        case 'ode'
            [pit,t] = ctmc_transient(InfGen,pi0,options.timespan(1),options.timespan(2),options.stiff,[],options.timestep);
        otherwise
            line_error(mfilename,sprintf('Unknown options.config.transient_method ''%s''; use ''ode'' or ''fau''.',transient_method));
    end
    mscale = [];
else
    if ~strcmp(transient_method,'ode')
        % A time-INHOMOGENEOUS generator is integrated by its own
        % piecewise-constant propagator below; accepting 'fau' here would
        % silently answer the constant-rate question instead.
        line_error(mfilename,'options.config.transient_method is not available together with a rate schedule.');
    end
    [pit,t,mscale] = local_ctmc_timevarying(sn, options, InfGen, StateSpace, pi0, rate_sched, M, K);
end
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
        % see _kb/06-solver-catalog.md (CTMC section, time-inhomogeneous generator) for rationale
        if ~isempty(mscale)
            TNt{ist,k} = TNt{ist,k} .* squeeze(mscale(ist,k,:));
        end
        % see _kb/06-solver-catalog.md (CTMC section, time-inhomogeneous generator) for rationale
        if sn.nodetype(sn.stationToNode(ist)) == NodeType.Source
            QNt{ist,k} = zeros(size(pit,1),1);
            UNt{ist,k} = zeros(size(pit,1),1);
            continue
        end
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

function [pit, t] = local_ctmc_fau(Q, pi0, options)
% Transient trajectory by FAST ADAPTIVE UNIFORMIZATION, selected with
% options.config.transient_method = 'fau' (see CTMC_FAU for the method).
%
% WHY IT IS NOT A METHOD NAME. 'fau' is not in SolverCTMC.listValidMethods and
% must not be: that list is enumerated by the sanity harness, which then demands
% a recorded baseline per method, and this one changes no stationary answer at
% all -- it is the transient path only. The config key is where the other
% transient switches of this analyzer already live (rate_sched, ctmc_tv_ngrid).
%
% MARCHED, NOT RESTARTED. pi(t_{k+1}) is obtained from pi(t_k) over the step
% rather than from pi(0) over the whole horizon, which is what keeps the cost
% proportional to the grid instead of quadratic in it. Every step removes a
% little mass and none puts any back, so the per-step tolerance is EPSILON
% divided by the number of steps and the total defect stays below EPSILON; the
% accumulated defect is reported rather than normalized away, since the whole
% point of the method is that its error is a measured quantity.
ts = options.timespan;
if ~isfinite(ts(2))
    line_error(mfilename,'transient_method ''fau'' needs a finite horizon; options.timespan(2) is not finite.');
end
if isfield(options,'timestep') && ~isempty(options.timestep) && options.timestep > 0
    t = (ts(1):options.timestep:ts(2))';
    if t(end) ~= ts(2)
        t = [t; ts(2)];
    end
else
    ngrid = 100;
    if isfield(options,'config') && isfield(options.config,'fau_ngrid') && ~isempty(options.config.fau_ngrid)
        ngrid = options.config.fau_ngrid;
    end
    t = linspace(ts(1), ts(2), ngrid)';
end
nt = numel(t);
epsilon = 1e-6;
if isfield(options,'config') && isfield(options.config,'fau_epsilon') && ~isempty(options.config.fau_epsilon)
    epsilon = options.config.fau_epsilon;
end
delta = 1e-12;
if isfield(options,'config') && isfield(options.config,'fau_delta') && ~isempty(options.config.fau_delta)
    delta = options.config.fau_delta;
end
epsStep = epsilon / max(1, nt-1);

pit = zeros(nt, length(Q));
pit(1,:) = pi0(:)';
defect = 0;
suppmax = 0;
for k = 1:nt-1
    [pk, info] = ctmc_fau(pit(k,:), Q, t(k+1)-t(k), epsStep, delta);
    pit(k+1,:) = pk;
    defect = defect + info.errorBound;
    suppmax = max(suppmax, info.supportMax);
end
line_debug('CTMC transient by FAU: %d grid points, support at most %d of %d states, missing mass %.3e', ...
    nt, suppmax, length(Q), defect);
if defect > GlobalConstants.CoarseTol
    line_warning(mfilename,'FAU transient discarded %.3e of the probability mass over the horizon; tighten options.config.fau_epsilon or options.config.fau_delta.\n', defect);
end
end

function [pit, t, mscale] = local_ctmc_timevarying(sn, options, Qbase, StateSpace, pi0, rate_sched, M, K)
% Integrate the time-inhomogeneous forward equation dpi/dt = pi Q(t) by
% non-homogeneous uniformization over a piecewise-constant grid. The generator
% is Q(t) = Qbase + sum_sc (m_sc(t)-1) Qhat_sc, where Qhat_sc is the linear
% component of Qbase attributable to the scaled (station,class) rate,
% extracted by a single probe rebuild (Q is linear in sn.rates).
ts = options.timespan;
Ngrid = 100;
if isfield(options,'config') && isfield(options.config,'ctmc_tv_ngrid') && ~isempty(options.config.ctmc_tv_ngrid)
    Ngrid = options.config.ctmc_tv_ngrid;
end
t = linspace(ts(1), ts(2), Ngrid)';
nt = numel(t);
nS = size(Qbase,1);

% Build the per-schedule generator component and multiplier trajectory.
probe = 2.0;
nsc = numel(rate_sched);
Qhat = cell(1,nsc);
mtraj = ones(nt, nsc);
scStation = zeros(1,nsc); scClass = zeros(1,nsc);
for s = 1:nsc
    ist = rate_sched(s).station;
    r = rate_sched(s).class;
    scStation(s) = ist; scClass(s) = r;
    % see _kb/06-solver-catalog.md (CTMC section, time-inhomogeneous generator) for rationale
    snp = sn;
    snp.rates(ist,r) = snp.rates(ist,r) * probe;
    if iscell(snp.proc) && numel(snp.proc) >= ist && iscell(snp.proc{ist}) && numel(snp.proc{ist}) >= r
        pr = snp.proc{ist}{r};
        if iscell(pr)
            for z = 1:numel(pr)
                if isnumeric(pr{z})
                    pr{z} = pr{z} * probe;
                end
            end
            snp.proc{ist}{r} = pr;
        end
    end
    if iscell(snp.mu) && numel(snp.mu) >= ist && iscell(snp.mu{ist}) && numel(snp.mu{ist}) >= r
        snp.mu{ist}{r} = snp.mu{ist}{r} * probe;
    end
    Qp = solver_ctmc(snp, options);
    if size(Qp,1) ~= nS
        line_error(mfilename, 'rate_sched probe changed the CTMC state-space size; cannot build time-varying generator.');
    end
    Qhat{s} = (Qp - Qbase) / (probe - 1);
    % multiplier m(t) = rate(t)/nominal (nominal defaults to sn.rates(ist,r))
    if isfield(rate_sched(s),'nominal') && ~isempty(rate_sched(s).nominal)
        nominal = rate_sched(s).nominal;
    else
        nominal = sn.rates(ist,r);
    end
    seg_t = rate_sched(s).tgrid(:)';
    seg_r = rate_sched(s).rates(:)';
    mtraj(:,s) = interp1(seg_t, seg_r, min(max(t,seg_t(1)),seg_t(end)), 'linear') / nominal;
end

% see _kb/06-solver-catalog.md (CTMC section, time-inhomogeneous generator) for rationale
pit = zeros(nt, nS);
pit(1,:) = pi0(:)';
Qfull = full(Qbase);
QhatFull = cell(1,nsc);
for s = 1:nsc
    QhatFull{s} = full(Qhat{s});
end
for k = 1:nt-1
    dt = t(k+1) - t(k);
    Qk = Qfull;
    for s = 1:nsc
        mk = 0.5*(mtraj(k,s) + mtraj(k+1,s));
        Qk = Qk + (mk - 1) * QhatFull{s};
    end
    pit(k+1,:) = pit(k,:) * expm(Qk * dt);
end

% Per-(station,class) throughput multiplier over time (1 where not scaled).
mscale = ones(M, K, nt);
for s = 1:nsc
    mscale(scStation(s), scClass(s), :) = reshape(mtraj(:,s), 1, 1, nt);
end
end
