function [tgrid, Mmat, tbreaks] = solver_fluid_ratemult(numEvents, M, K, enabled, q_indices, Kic, Mu, eventIdx, options)
% [TGRID, MMAT, TBREAKS] = SOLVER_FLUID_RATEMULT(NUMEVENTS, M, K, ENABLED, Q_INDICES, KIC, MU, EVENTIDX, OPTIONS)
%
% Build the time-varying per-event rate multiplier for the closing fluid ODE.
% Returns TGRID (1 x ngrid, strictly increasing) and MMAT (numEvents x ngrid);
% event e is scaled at time t by fluid_interpcols(TGRID, MMAT, t)(e). Returns
% [] when no time-varying source is configured, so the caller keeps the legacy
% autonomous closure.
%
% TBREAKS lists the instants at which the multiplier JUMPS, i.e. the NHPP
% schedule's own segment bounds. They are reported separately from TGRID
% because they are the only grid instants an integrator must not step across:
% the schedule is piecewise CONSTANT, so the drift is discontinuous there and
% nothing in the right-hand side tells a step controller where the jump is.
% SOLVER_FLUID_ITERATION makes each one an integration boundary. The other two
% sources carry no breaks: RATE_TRAJ and RATE_SCHED are sampled trajectories
% meant to be read as piecewise linear, which every integrator handles.
%
% Two independent, composable sources are honoured (both reduce to a per-event
% multiplicative factor because the closing rate is rate = rateBase .* theta(x)
% and rateBase is linear in the station-class service/arrival rate):
%   (1) options.config.rate_traj = {tgrid, Mmat} - a caller-supplied event
%       multiplier matrix (used by the coupled LN layer transient).
%   (2) options.config.nhpp_sched - a struct array of non-homogeneous
%       (NHPP) source intensities, each with fields:
%         .station  station index in the sn station space (must be EXT/source)
%         .class    class index
%         .nhpp     the process handle exposing getRateAt(t) and, for grid
%                   construction, getBreakpoints/getPeriod/isCyclic
%       The nominal (time-average) rate baked into rateBase for that
%       station-class is Mu{station}{class}(1); the multiplier is
%       getRateAt(t)/nominal, applied to every event sourced at that station
%       and class (eventIdx == q_indices(station,class) + kic - 1).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tgrid = [];
Mmat = [];
tbreaks = [];
if ~isfield(options,'config') || isempty(options.config)
    return
end
cfg = options.config;

% -- source (1): caller-supplied rate_traj ------------------------------------
user_tg = [];
user_M = [];
if isfield(cfg,'rate_traj') && ~isempty(cfg.rate_traj)
    user_tg = cfg.rate_traj{1}(:)';
    user_M = cfg.rate_traj{2};
    if size(user_M,1) ~= numEvents
        line_error(mfilename, sprintf('rate_traj multiplier matrix has %d rows but the closing ODE has %d events.', size(user_M,1), numEvents));
    end
end

% -- source (2): NHPP source intensities --------------------------------------
nhpp_tg = [];
nhpp_M = [];
if isfield(cfg,'nhpp_sched') && ~isempty(cfg.nhpp_sched)
    sched = cfg.nhpp_sched;
    % horizon over which to expand a (possibly cyclic) schedule
    t0 = 0;
    tend = Inf;
    if isfield(options,'timespan') && numel(options.timespan) >= 2
        if isfinite(options.timespan(1)), t0 = options.timespan(1); end
        tend = options.timespan(2);
    end
    for s = 1:numel(sched)
        i = sched(s).station;
        c = sched(s).class;
        nh = sched(s).nhpp;
        if ~enabled(i,c)
            continue % class not served/active at this station in the fluid ODE
        end
        nominal = Mu{i}{c}(1);
        if ~(nominal > 0)
            continue
        end
        % horizon: for a non-finite timespan use a few periods so a cyclic
        % schedule is represented rather than clamped after one segment
        period = nh.getPeriod();
        if ~isfinite(tend)
            if isfinite(period) && period > 0
                thi = t0 + 3*period;
            else
                thi = t0 + 1;
            end
        else
            thi = tend;
        end
        [seg_t, seg_r, seg_b] = local_nhpp_steps(nh, t0, thi);
        tbreaks = [tbreaks, seg_b]; %#ok<AGROW>
        rowmult = seg_r / nominal;
        % rows: all events sourced at (i,c) across its phases
        rows = false(numEvents,1);
        for kic = 1:Kic(i,c)
            rows = rows | (eventIdx(:) == (q_indices(i,c) + kic - 1));
        end
        thisM = ones(numEvents, numel(seg_t));
        thisM(rows,:) = repmat(rowmult, sum(rows), 1);
        [nhpp_tg, nhpp_M] = local_merge(nhpp_tg, nhpp_M, seg_t, thisM, numEvents);
    end
end

% -- source (3): explicit per-(station,class) rate trajectories ---------------
% options.config.rate_sched is a struct array with fields:
%   .station, .class  station/class in the sn space
%   .tgrid, .rates    the time-varying rate (same length)
%   .nominal (opt)    nominal rate baked into rateBase (default Mu{i}{c}(1))
% Used by the coupled LN layer transient to inject time-varying inter-layer
% demand (think-time / service) trajectories, reusing the same station-class
% -> event expansion as the NHPP path.
sched_tg = [];
sched_M = [];
if isfield(cfg,'rate_sched') && ~isempty(cfg.rate_sched)
    rs = cfg.rate_sched;
    for s = 1:numel(rs)
        i = rs(s).station;
        c = rs(s).class;
        if ~enabled(i,c)
            continue
        end
        if isfield(rs(s),'nominal') && ~isempty(rs(s).nominal)
            nominal = rs(s).nominal;
        else
            nominal = Mu{i}{c}(1);
        end
        if ~(nominal > 0)
            continue
        end
        seg_t = rs(s).tgrid(:)';
        rowmult = rs(s).rates(:)' / nominal;
        rows = false(numEvents,1);
        for kic = 1:Kic(i,c)
            rows = rows | (eventIdx(:) == (q_indices(i,c) + kic - 1));
        end
        thisM = ones(numEvents, numel(seg_t));
        thisM(rows,:) = repmat(rowmult, sum(rows), 1);
        [sched_tg, sched_M] = local_merge(sched_tg, sched_M, seg_t, thisM, numEvents);
    end
end

% -- compose all sources ------------------------------------------------------
[tgrid, Mmat] = local_merge(user_tg, user_M, nhpp_tg, nhpp_M, numEvents);
[tgrid, Mmat] = local_merge(tgrid, Mmat, sched_tg, sched_M, numEvents);
tbreaks = unique(tbreaks(:)');
end

function [seg_t, seg_r, seg_b] = local_nhpp_steps(nh, t0, thi)
% Build a step-faithful (time, rate) sampling of a piecewise-constant NHPP
% intensity over [t0, thi]. Each segment contributes two samples at its start
% and just before its end, so clamped-linear interpolation reproduces the step
% with a negligible transition ramp. SEG_B returns the INTERIOR segment bounds,
% i.e. where that ramp sits and the intensity actually jumps.
bp = nh.getBreakpoints();
bp = bp(:)';
period = nh.getPeriod();
if nh.isCyclic() && isfinite(period) && period > 0
    bounds = [];
    kmax = ceil((thi - t0)/period) + 2;
    for k = -1:kmax
        bounds = [bounds, bp + k*period]; %#ok<AGROW>
    end
else
    bounds = bp;
end
bounds = unique([t0, bounds(bounds > t0 & bounds < thi), thi]);
neps = max(1e-9, 1e-6*(thi - t0));
nb = numel(bounds) - 1;
seg_t = zeros(1, 2*nb);
seg_r = zeros(1, 2*nb);
for k = 1:nb
    a = bounds(k);
    b = bounds(k+1);
    r = nh.getRateAt((a + b)/2);
    seg_t(2*k-1) = a;
    seg_r(2*k-1) = r;
    seg_t(2*k) = max(a + neps, b - neps);
    seg_r(2*k) = r;
end
seg_b = bounds(2:end-1);
end

function [tg, Mg] = local_merge(tg1, M1, tg2, M2, numEvents)
% Merge two event-multiplier trajectories onto the union time grid by
% elementwise product (identity where a source is silent).
if isempty(M1)
    tg = tg2; Mg = M2; return
end
if isempty(M2)
    tg = tg1; Mg = M1; return
end
tg = unique([tg1(:)', tg2(:)']);
ng = numel(tg);
Mg = ones(numEvents, ng);
for j = 1:ng
    Mg(:,j) = fluid_interpcols(tg1, M1, tg(j)) .* fluid_interpcols(tg2, M2, tg(j));
end
end
