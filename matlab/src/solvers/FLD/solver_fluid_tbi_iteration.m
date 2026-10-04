function [xvec_it, xvec_t, t, iter] = solver_fluid_tbi_iteration(sn, N, Mu, Phi, PH, P, S, xvec_it, ydefault, slowrate, Tstart, max_time, options)
% [XVEC_IT, XVEC_T, T, ITER] = SOLVER_FLUID_TBI_ITERATION(SN, N, MU, PHI, PH, P, S, XVEC_IT, YDEFAULT, SLOWRATE, TSTART, MAX_TIME, OPTIONS)
%
% Trajectory-based iteration (TBI) for the transient fluid solution.
% The station set is partitioned into cells (tbi_partition). On each time
% segment, the IVP of every cell is solved with the state of the other
% cells frozen at the trajectory computed in the previous iteration (Gauss-Seidel or Jacobi
% waveform relaxation); iterations repeat until the trajectory sup-norm gap
% falls below options.config.tbi_tol. Cross-cell inflows are therefore
% evaluated on frozen trajectories, interior flows on the live cell state,
% matching the decomposed ODEs of the TBI method.
%
% Reference: Sheldon, Tuncer, Casale, "TBI: Transient Hierarchical
% Modeling of Large-Scale Vehicle Sharing Systems", IEEE T-ITS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter_max = options.iter_max;
tol = options.tol;
timespan = options.timespan;
stiff = options.stiff;

% TBI configuration
if isfield(options.config,'tbi_tol') && ~isempty(options.config.tbi_tol)
    tbi_tol = options.config.tbi_tol;
else
    tbi_tol = 1e-3; % absolute sup-norm tolerance on state trajectories
end
if isfield(options.config,'tbi_iter_max') && ~isempty(options.config.tbi_iter_max)
    tbi_iter_max = options.config.tbi_iter_max;
else
    tbi_iter_max = 50; % maximum TBI iterations per time segment
end
if isfield(options.config,'hide_immediate') && options.config.hide_immediate
    line_warning(mfilename,'The tbi method does not support hide_immediate, ignoring it.\n');
end
if isfield(options.config,'tbi_parallel') && ~isempty(options.config.tbi_parallel)
    tbi_parallel = options.config.tbi_parallel && license('test','Distrib_Computing_Toolbox');
    if options.config.tbi_parallel && ~tbi_parallel
        line_warning(mfilename,'tbi_parallel requested but Parallel Computing Toolbox is unavailable, running serially.\n');
    end
else
    tbi_parallel = false;
end
% iteration ordering: Gauss-Seidel (paper Alg. 1, cells see the trajectories of
% cells already solved in the same iteration) or Jacobi (all cells see the
% previous iteration; required by the parallel path)
if isfield(options.config,'tbi_iteration') && ~isempty(options.config.tbi_iteration)
    tbi_gs = strcmpi(options.config.tbi_iteration,'gs') || strcmpi(options.config.tbi_iteration,'gauss-seidel');
else
    tbi_gs = true; % default
end
if tbi_parallel && tbi_gs
    tbi_gs = false; % Gauss-Seidel is sequential by construction
end

M = sn.nstations;
K = sn.nclasses;
sched = sn.sched;
schedparam = sn.schedparam;

% state indexing, mirrors solver_fluid_odes
w = ones(M,K);
enabled = false(M,K);
q_indices = zeros(M,K);
Kic = zeros(M,K);
csum = 1;
for i = 1:M
    for c = 1:K
        if isnan(Mu{i}{c})
            numphases = 0;
            q_indices(i,c) = csum;
        elseif isempty(Mu{i}{c})
            numphases = 0;
            q_indices(i,c) = csum;
        else
            numphases = length(Mu{i}{c});
            q_indices(i,c) = csum;
            enabled(i,c) = true;
        end
        Kic(i,c) = numphases;
        csum = csum + numphases;
    end
    if sched(i) == SchedStrategy.DPS
        w(i,:) = schedparam(i,:);
    end
end
ndim = csum - 1;

% closing-method event machinery, built once on the full model
all_jumps = ode_jumps_new(M, K, enabled, q_indices, P, Kic);
[rateBase, eventIdx] = ode_rate_base(sn, Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, sched, all_jumps);
% Same reduction the other routes take: a cell of the trajectory-based
% partition integrates the same drift, so an InfRate coordinate costs it the
% same and is complemented away here too.
if fluid_hide_immediate(sn, options)
    [all_jumps, rateBase, eventIdx] = ...
        ode_eliminate_immediate(all_jumps, rateBase, eventIdx, sn, options);
end
rates_h = @(x) ode_rates_closing(x, M, K, enabled, q_indices, Kic, S, w, sched, rateBase, eventIdx);

% partition stations into cells and build per-cell state masks; the
% restriction of the jump matrix to a cell's rows retains interior and
% boundary events only through their effect on the cell state, so inbound
% events read the frozen complement while outbound mass leaves the cell
cells = tbi_partition(sn, options);
ncells = numel(cells);
cellmask = cell(1,ncells);
Jint = cell(1,ncells);   % interior/outbound jumps restricted to cell rows
Jext = cell(1,ncells);   % inbound jumps sourced at complement stations
Eext = cell(1,ncells);   % events feeding Jext (rates frozen per iteration)
Eint = cell(1,ncells);   % events sourced inside the cell (global indices)
Mk = zeros(1,ncells);
enabled_k = cell(1,ncells); q_idx_k = cell(1,ncells); Kic_k = cell(1,ncells);
S_k = cell(1,ncells); w_k = cell(1,ncells); sched_k = cell(1,ncells);
rateBase_k = cell(1,ncells); eventIdx_k = cell(1,ncells);
for kc = 1:ncells
    cs = sort(cells{kc}(:)');
    % local station-restricted structures; the local state layout follows
    % ascending global index order, so it coincides with the sorted mask
    Mk(kc) = numel(cs);
    enabled_k{kc} = enabled(cs,:);
    Kic_k{kc} = Kic(cs,:);
    S_k{kc} = S(cs);
    w_k{kc} = w(cs,:);
    sched_k{kc} = sched(cs);
    q_idx_k{kc} = zeros(Mk(kc),K);
    lsum = 1;
    for a = 1:Mk(kc)
        for c = 1:K
            q_idx_k{kc}(a,c) = lsum;
            lsum = lsum + Kic_k{kc}(a,c);
        end
    end
    mask = false(ndim,1);
    for i = cs
        for c = 1:K
            if Kic(i,c) > 0
                mask(q_indices(i,c):(q_indices(i,c)+Kic(i,c)-1)) = true;
            end
        end
    end
    idx = find(mask);
    cellmask{kc} = idx;
    glob2loc = zeros(ndim,1);
    glob2loc(idx) = 1:numel(idx);
    % events sourced inside the cell: rates depend only on the live cell
    % state because the closing scaling is local to the source station
    isInt = glob2loc(eventIdx) > 0;
    Eint{kc} = find(isInt);
    rateBase_k{kc} = rateBase(Eint{kc});
    eventIdx_k{kc} = glob2loc(eventIdx(Eint{kc}));
    Jint{kc} = all_jumps(idx, Eint{kc});
    % inbound events: sourced at complement stations but move mass into the
    % cell; their rates are frozen per iteration and precomputed on the grid
    Ee = find(~isInt);
    Je = all_jumps(idx, Ee);
    keep = any(Je ~= 0, 1);
    Eext{kc} = Ee(keep);
    Jext{kc} = Je(:, keep);
end

% horizon growth heuristic, mirrors solver_fluid_iteration
nonZeroRates = slowrate(:);
nonZeroRates = nonZeroRates(nonZeroRates > tol);
nonZeroRates = nonZeroRates(isfinite(nonZeroRates));
if isempty(nonZeroRates)
    nonZeroRates = 1; % fallback when all rates are zero or infinite
end

goon = true;
iter = 0;
t = [];
xvec_t = [];
T0 = timespan(1);
T = 0;
while (isfinite(timespan(2)) && T < timespan(2)) || (goon && iter < iter_max)
    iter = iter + 1;
    if toc(Tstart) > max_time
        goon = false;
        break;
    end

    y0 = xvec_it{iter-1 +1}(:)';
    if iter == 1 % first iteration
        T = min(timespan(2),abs(10/min(nonZeroRates)));
    else
        T = min(timespan(2),abs(10*iter/min(nonZeroRates)));
    end
    trange = [T0, T];

    % frozen trajectory on this segment, initialized constant at the
    % segment entry state (warm start of the waveform relaxation)
    tprev = [T0; T];
    Yprev = [y0; y0];

    delta = Inf;
    for itn = 1:tbi_iter_max
        % precompute the full closing-rate vector on the frozen grid once
        % per iteration (shared by all cells), then collapse the inbound events
        % of each cell to a drift time series of cell dimension
        ngrid = numel(tprev);
        Rfull = zeros(numel(rateBase), ngrid);
        for jg = 1:ngrid
            Rfull(:,jg) = rates_h(Yprev(jg,:)');
        end
        cell_t = cell(1,ncells);
        cell_y = cell(1,ncells);
        if tbi_parallel
            parfor kc = 1:ncells
                [cell_t{kc}, cell_y{kc}] = tbi_solve_cell(kc, tprev, Yprev, Rfull, trange, y0, ydefault, ...
                    cellmask, Jint, Jext, Eext, Mk, K, enabled_k, q_idx_k, Kic_k, S_k, w_k, sched_k, ...
                    rateBase_k, eventIdx_k, tol, stiff, options);
            end
        else
            for kc = 1:ncells
                [cell_t{kc}, cell_y{kc}] = tbi_solve_cell(kc, tprev, Yprev, Rfull, trange, y0, ydefault, ...
                    cellmask, Jint, Jext, Eext, Mk, K, enabled_k, q_idx_k, Kic_k, S_k, w_k, sched_k, ...
                    rateBase_k, eventIdx_k, tol, stiff, options);
                if tbi_gs
                    % Gauss-Seidel: refresh this cell's event rates on the
                    % exchange grid so later cells see the new trajectory
                    tc = cell_t{kc};
                    yl = interp1(tc, cell_y{kc}, min(max(tprev, tc(1)), tc(end)));
                    for jg = 1:ngrid
                        Rfull(Eint{kc}, jg) = ode_rates_closing(yl(jg,:)', Mk(kc), K, enabled_k{kc}, q_idx_k{kc}, Kic_k{kc}, S_k{kc}, w_k{kc}, sched_k{kc}, rateBase_k{kc}, eventIdx_k{kc});
                    end
                end
            end
        end
        tgrid = tprev;
        for kc = 1:ncells
            tgrid = union(tgrid, cell_t{kc});
        end
        % assemble the new full trajectory on the union time grid
        Ynew = zeros(numel(tgrid), ndim);
        for kc = 1:ncells
            tc = cell_t{kc};
            tq = min(max(tgrid, tc(1)), tc(end));
            Ynew(:,cellmask{kc}) = interp1(tc, cell_y{kc}, tq);
        end
        % sup-norm gap against the previous iteration trajectory
        tq = min(max(tgrid, tprev(1)), tprev(end));
        Yold = interp1(tprev, Yprev, tq);
        delta = max(abs(Ynew(:) - Yold(:)));
        tprev = tgrid;
        Yprev = Ynew;
        if delta < tbi_tol || toc(Tstart) > max_time
            break
        end
    end
    if delta >= tbi_tol && options.verbose > 0
        line_warning(mfilename,'TBI iterations did not converge within tbi_iter_max=%d on segment [%g,%g], residual gap %g.\n', tbi_iter_max, T0, T, delta);
    end
    if options.verbose >= 2
        line_printf('\nTBI segment %d [%g,%g]: %d iterations, %d cell solves, grid %d points, gap %g', iter, T0, T, itn, itn*ncells, numel(tprev), delta);
    end

    xvec_t(end+1:end+size(Yprev,1),:) = Yprev;
    t(end+1:end+numel(tprev),:) = tprev;
    xvec_it{iter +1} = xvec_t(end,:);
    if options.verbose >= 2 && isfinite(N)
        % closed-population conservation check; waveform relaxation loses
        % mass transiently, at convergence the total must return to N
        massgap = abs(sum(xvec_it{iter +1}) - N);
        if massgap > max(0.01*N, 10*tbi_tol)
            line_warning(mfilename,'TBI mass conservation gap %g at t=%g (N=%g).\n', massgap, T, N);
        end
    end
    T0 = T; % for next segment

    if T >= timespan(2)
        goon = false;
    end
end
end

function [tc, yc] = tbi_solve_cell(kc, tprev, Yprev, Rfull, trange, y0, ydefault, ...
    cellmask, Jint, Jext, Eext, Mk, K, enabled_k, q_idx_k, Kic_k, S_k, w_k, sched_k, ...
    rateBase_k, eventIdx_k, tol, stiff, options)
% solve one cell IVP against the frozen complement (one relaxation step)
idx = cellmask{kc};
Bk = Jext{kc} * Rfull(Eext{kc},:); % inbound drift on the frozen grid
ode_c = @(tt,xc) Jint{kc} * ode_rates_closing(xc, Mk(kc), K, enabled_k{kc}, q_idx_k{kc}, Kic_k{kc}, S_k{kc}, w_k{kc}, sched_k{kc}, rateBase_k{kc}, eventIdx_k{kc}) ...
    + fluid_interpcols(tprev, Bk, tt);
odeopt_c = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:numel(idx));
try
    if stiff
        [tc, yc] = ode_solve_stiff(ode_c, trange, y0(idx), odeopt_c, options);
    else
        [tc, yc] = ode_solve(ode_c, trange, y0(idx), odeopt_c, options);
    end
catch
    line_printf('\nThe initial point is invalid, Fluid solver switching to default initialization.');
    [tc, yc] = ode_solve(ode_c, trange, ydefault(idx), odeopt_c, options);
end
end

% Cross-cell inflow interpolation is provided by the shared helper
% fluid_interpcols (formerly the local tbi_interpcols).
