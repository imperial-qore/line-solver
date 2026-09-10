function [QNlqn_t, UNlqn_t, TNlqn_t] = getTranAvgCoupled(self, Qt, Ut, Tt) %#ok<INUSD>
% [QNLQN_T,UNLQN_T,TNLQN_T] = GETTRANAVGCOUPLED(SELF,QT,UT,TT)
%
% Coupled layered transient by waveform relaxation over the LQN ensemble.
% Unlike getTranAvgDecoupled, which freezes inter-layer demands at the
% converged fixed point, this reconciles the per-layer transients iteratively:
% each layer's fluid transient is driven by TIME-VARYING inter-layer demand
% trajectories taken from the other layers' latest transients, and the loop
% repeats until the trajectories stop changing (sup-norm gap over time). The
% time-varying demands are injected into each layer's closing ODE through the
% per-event rate multiplier (options.config.rate_sched -> solver_fluid_ratemult).
%
% Iteration 0 uses the frozen equilibrium demands, so it reproduces
% getTranAvgDecoupled exactly; at convergence every layer relaxes to its fixed
% point, so the endpoint equals getAvg. The return layout is the same
% block-diagonal (station x class per layer) as getTranAvgDecoupled.
%
% Coupled channels: task think times (client delay) and synchronous-call
% service demands (caller client station). Both are the dominant inter-layer
% couplings; intra-layer host service stays at its equilibrium value.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lqn = self.lqn;
% Capture the transient horizon BEFORE getAvg.
% see _kb/06-solver-catalog.md (LN section) for rationale
ts = self.options.timespan;
self.getAvg; % converge the fixed point; layer solvers primed with equilibrium

if ~(numel(ts) >= 2 && all(isfinite(ts)))
    % No finite transient horizon: nothing to co-evolve, defer to decoupled.
    [QNlqn_t, UNlqn_t, TNlqn_t] = self.getTranAvgDecoupled();
    return
end

E = self.nlayers;

% relaxation controls (reuse the ensemble iteration budget / tolerance)
maxit = 20;
if isfield(self.options,'config') && isfield(self.options.config,'ln_transient_iter_max') ...
        && ~isempty(self.options.config.ln_transient_iter_max)
    maxit = self.options.config.ln_transient_iter_max;
end
tol = 1e-2;
if isfield(self.options,'config') && isfield(self.options.config,'ln_transient_tol') ...
        && ~isempty(self.options.config.ln_transient_tol)
    tol = self.options.config.ln_transient_tol;
end
Ngrid = 100;
tgrid = linspace(ts(1), ts(2), Ngrid)';

% Per-layer sn (for node->station and class bookkeeping), cached once.
layerSn = cell(1,E);
for e = 1:E
    layerSn{e} = self.ensemble{e}.getStruct;
end

% Iteration 0: decoupled transients (frozen equilibrium demands already set by
% getAvg). local_run_layers returns per-layer {M,K} structs on the layer's own
% time grid, plus per-(station,class) trajectories resampled onto tgrid.
[blocks, traj] = local_run_layers(self, ts, tgrid, cell(1,E), layerSn);

for iter = 1:maxit
    trajPrev = traj;
    % 1) recompute inter-layer demand trajectories from the latest layer traj
    demand = local_recompute_demand(self, traj, tgrid, layerSn);
    % 2) build the per-layer rate_sched injections from those demands
    schedByLayer = local_build_rate_sched(self, demand, tgrid, layerSn);
    % 3) re-run each layer with its injected time-varying demand
    [blocks, traj] = local_run_layers(self, ts, tgrid, schedByLayer, layerSn);
    % 4) convergence: sup-norm gap of the queue-length trajectories
    gap = local_supnorm_gap(trajPrev, traj);
    if self.options.verbose >= VerboseLevel.STD
        line_printf('\nLN coupled transient: iter %d, sup-norm gap %.3e', iter, gap);
    end
    if gap < tol
        break
    end
end

% Assemble the block-diagonal aggregate exactly as getTranAvgDecoupled does.
QNlqn_t = cell(0,0);
UNlqn_t = cell(0,0);
TNlqn_t = cell(0,0);
for e = 1:E
    [crows, ccols] = size(QNlqn_t);
    Qe = blocks{e}.Q; Ue = blocks{e}.U; Te = blocks{e}.T;
    QNlqn_t(crows+1:crows+size(Qe,1), ccols+1:ccols+size(Qe,2)) = Qe;
    UNlqn_t(crows+1:crows+size(Ue,1), ccols+1:ccols+size(Ue,2)) = Ue;
    TNlqn_t(crows+1:crows+size(Te,1), ccols+1:ccols+size(Te,2)) = Te;
end
end

% -------------------------------------------------------------------------
function [blocks, traj] = local_run_layers(self, ts, tgrid, schedByLayer, layerSn)
% Run each layer's transient (optionally with an injected rate_sched) and
% return: blocks{e} = struct with {M,K} cell fields Q,U,T (native handles, for
% block-diagonal assembly); traj{e} = struct with M x K x numel(tgrid) arrays
% Q,U,T,R resampled onto tgrid (R = Q./T residence via Little).
E = self.nlayers;
blocks = cell(1,E);
traj = cell(1,E);
ng = numel(tgrid);
for e = 1:E
    s = self.solvers{e};
    savedTs = s.options.timespan;
    hadSched = isfield(s.options.config,'rate_sched');
    if hadSched
        savedSched = s.options.config.rate_sched;
    end
    s.options.timespan = ts;
    if ~isempty(schedByLayer{e})
        s.options.config.rate_sched = schedByLayer{e};
    elseif hadSched
        s.options.config.rate_sched = [];
    end
    % Drop any cached transient before re-solving with injected rate_sched.
    % see _kb/06-solver-catalog.md (LN section) for rationale
    s.reset();
    [Qe, Ue, Te] = s.getTranAvg();
    s.options.timespan = savedTs;
    if hadSched
        s.options.config.rate_sched = savedSched;
    elseif isfield(s.options.config,'rate_sched')
        s.options.config = rmfield(s.options.config,'rate_sched');
    end
    be = struct();       % avoid struct('Q',{cell}) which builds a struct ARRAY
    be.Q = Qe; be.U = Ue; be.T = Te;
    blocks{e} = be;
    M = size(Qe,1); K = size(Qe,2);
    Q = zeros(M,K,ng); U = zeros(M,K,ng); T = zeros(M,K,ng); R = zeros(M,K,ng);
    for i = 1:M
        for r = 1:K
            qv = local_resample(Qe{i,r}, tgrid);
            uv = local_resample(Ue{i,r}, tgrid);
            tv = local_resample(Te{i,r}, tgrid);
            Q(i,r,:) = qv; U(i,r,:) = uv; T(i,r,:) = tv;
            R(i,r,:) = qv ./ max(tv, GlobalConstants.FineTol);
        end
    end
    traj{e} = struct('Q',Q,'U',U,'T',T,'R',R);
end
end

% -------------------------------------------------------------------------
function v = local_resample(h, tgrid)
% Resample a transient handle (struct with .t/.metric, or a scalar) onto tgrid.
if isstruct(h) && isfield(h,'t') && isfield(h,'metric') && numel(h.t) >= 2
    v = interp1(h.t(:), h.metric(:), tgrid, 'linear', 'extrap');
elseif isstruct(h) && isfield(h,'metric')
    v = repmat(h.metric(end), numel(tgrid), 1);
else
    v = zeros(numel(tgrid),1);
end
end

% -------------------------------------------------------------------------
function demand = local_recompute_demand(self, traj, tgrid, layerSn)
% Recompute the time-varying inter-layer demands from the layer trajectories,
% pointwise in t, mirroring the scalar updateThinkTimes / updateMetricsDefault
% formulas. Returns struct with:
%   thinkt : dictionary tidx -> (ng x 1) think-time trajectory
%   callservt : dictionary cidx -> (ng x 1) call service-time trajectory
lqn = self.lqn;
ng = numel(tgrid);
thinkt = configureDictionary('double','cell');
callservt = configureDictionary('double','cell');

% Task think times: from the task's own server-layer utilization/throughput.
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if isnan(self.idxhash(tidx)) || lqn.isref(tidx)
        continue
    end
    [e, sIdx] = self.layerOf(tidx);
    K = layerSn{e}.nclasses;
    Uti = zeros(ng,1); Tti = zeros(ng,1);
    for r = 1:K
        Uti = Uti + squeeze(traj{e}.U(sIdx,r,:));
        Tti = Tti + squeeze(traj{e}.T(sIdx,r,:));
    end
    njobs = max(self.njobs(tidx,:));
    % same closure as updateThinkTimes, so the same gate: only a reference
    % task's think time is a per-request delay -- see lqn_ref_thinktime
    userthink = lqn_ref_thinktime(lqn, tidx);
    Tsafe = max(Tti, GlobalConstants.FineTol);
    if lqn.sched(tidx) == SchedStrategy.INF
        tk = (njobs - Uti) ./ Tsafe - userthink;
    else
        tk = njobs .* abs(1 - Uti) ./ Tsafe - userthink;
    end
    tk = max(GlobalConstants.Zero, tk) + userthink; % total mean incl. user think
    thinkt{tidx} = tk;
end

% Synchronous-call service demands: callee entry response time * call mean.
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) ~= CallType.SYNC
        continue
    end
    eidx = lqn.callpair(cidx,2);       % callee entry
    tidx = lqn.parent(eidx);           % callee task
    if isnan(self.idxhash(tidx))
        continue
    end
    [e, sIdx] = self.layerOf(tidx);
    % response time of the callee at its server, summed over the entry classes
    K = layerSn{e}.nclasses;
    Rc = zeros(ng,1);
    for r = 1:K
        typ = self.ensemble{e}.classes{r}.attribute(1);
        if typ == LayeredNetworkElement.ENTRY && self.ensemble{e}.classes{r}.attribute(2) == eidx
            Rc = Rc + squeeze(traj{e}.R(sIdx,r,:));
        end
    end
    if all(Rc == 0)
        % fall back to the entry's activities response time
        Rc = squeeze(sum(traj{e}.R(sIdx,:,:),2));
    end
    callmean = 1.0;
    if ~isempty(lqn.callproc{cidx}) && isa(lqn.callproc{cidx},'Distribution')
        callmean = lqn.callproc{cidx}.getMean;
    end
    callservt{cidx} = Rc * callmean;
end

demand = struct('thinkt', thinkt, 'callservt', callservt);
end

% -------------------------------------------------------------------------
function schedByLayer = local_build_rate_sched(self, demand, tgrid, layerSn)
% Map the recomputed demand trajectories to per-layer rate_sched injections,
% using the same update maps updateLayers uses to place setService calls.
% options.config.ln_transient_channels selects which inter-layer coupling
% channels are injected: 'both' (default), 'thinkt' (client-delay only), or
% 'callservt' (synchronous-call service only). Used to isolate each channel's
% contribution to the coupled transient.
E = self.nlayers;
schedByLayer = cell(1,E);
lqn = self.lqn;
channels = 'both';
if isfield(self.options,'config') && isfield(self.options.config,'ln_transient_channels') ...
        && ~isempty(self.options.config.ln_transient_channels)
    channels = lower(self.options.config.ln_transient_channels);
end

% think-time channel (client delay of caller tasks)
if any(strcmp(channels, {'both','thinkt'}))
    map = self.thinkt_classes_updmap;
    for r = 1:size(map,1)
        idx = map(r,1); aidx = map(r,2); nodeidx = map(r,3); classidx = map(r,4);
        e = self.idxhash(idx);
        if isnan(e) || nodeidx ~= self.ensemble{e}.attribute.clientIdx
            continue
        end
        if lqn.type(aidx) == LayeredNetworkElement.TASK && lqn.sched(aidx) ~= SchedStrategy.REF ...
                && isKey(demand.thinkt, aidx)
            d = demand.thinkt{aidx};
            schedByLayer{e} = local_add_sched(schedByLayer{e}, layerSn{e}, nodeidx, classidx, tgrid, d);
        end
    end
end

% call-service channel (client station of caller for each sync call)
if any(strcmp(channels, {'both','callservt'}))
    map = self.call_classes_updmap;
    for c = 1:size(map,1)
        idx = map(c,1); cidx = map(c,2); nodeidx = map(c,3); classidx = map(c,4);
        e = self.idxhash(idx);
        if isnan(e) || nodeidx ~= self.ensemble{e}.attribute.clientIdx
            continue
        end
        if isKey(demand.callservt, cidx)
            d = demand.callservt{cidx};
            schedByLayer{e} = local_add_sched(schedByLayer{e}, layerSn{e}, nodeidx, classidx, tgrid, d);
        end
    end
end
end

% -------------------------------------------------------------------------
function sched = local_add_sched(sched, sn, nodeidx, classidx, tgrid, demand)
% Append one rate_sched entry that MODULATES the layer's equilibrium rate by
% the ratio of the transient demand to its steady-state (end-of-horizon) value:
%   effective_rate(t) = nominal * demand(end)/demand(t).
% Passing rate = 1/demand(t) and nominal_denominator = 1/demand(end) makes the
% multiplier = demand(end)/demand(t), which is exactly 1 at the horizon end, so
% the layer relaxes to its unmodified fixed point (oracle 1: endpoint == getAvg)
% regardless of any small mismatch between the fluid residence Q/T and the
% scalar equilibrium demand. During the transient the ratio modulates the rate.
ist = sn.nodeToStation(nodeidx);
if isnan(ist)
    return
end
d = demand(:)';
dend = d(end);
if ~(dend > GlobalConstants.FineTol)
    return % degenerate steady-state demand; skip this channel
end
% Bound the transient demand to a physical band around its steady-state value.
% see _kb/06-solver-catalog.md (LN section) for rationale
Cap = 20;
d = min(max(d, dend/Cap), dend*Cap);
% self-normalising rate: nominal is the reciprocal of the steady-state demand,
% so solver_fluid_ratemult's multiplier = rate/nominal = dend/d(t).
rates = 1 ./ d;
nominal = 1 / dend;
entry = struct('station', ist, 'class', classidx, 'tgrid', tgrid(:)', ...
    'rates', rates, 'nominal', nominal);
if isempty(sched)
    sched = entry;
else
    sched(end+1) = entry; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function gap = local_supnorm_gap(trajPrev, traj)
% Sup-norm of the queue-length trajectory change across all layers.
gap = 0;
for e = 1:numel(traj)
    d = abs(traj{e}.Q(:) - trajPrev{e}.Q(:));
    if ~isempty(d)
        gap = max(gap, max(d));
    end
end
end
