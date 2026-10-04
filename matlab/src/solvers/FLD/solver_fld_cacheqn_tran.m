function [tcache, hitprob_t, missprob_t, caches, arate, xocc, lamcell, sojourn] = solver_fld_cacheqn_tran(sn, options, x0cell, holdmap)
% SOLVER_FLD_CACHEQN_TRAN Transient refined-mean-field cache trajectory.
%
% Converges the per-class cache arrival rates with the same decomposition-
% aggregation alternation as the steady solver (da_cacheqn), then integrates
% the mean-field drift over OPTIONS.TIMESPAN to obtain the time-resolved
% per-class hit/miss probabilities of each cache. This is the transient
% counterpart of solver_fld_cacheqn_analyzer: the steady solver drives the
% drift to its fixed point, whereas here the same drift is integrated over
% the finite window from the supplied (or default) initial occupancy.
%
% SN:      NetworkStruct with at least one Cache node.
% OPTIONS: solver options; OPTIONS.TIMESPAN = [t0,t1] sets the window.
% X0CELL:  optional cell(1,ncaches); X0CELL{c} seeds cache c occupancy
%          (flat DDPP state vector, n_items*(h+1)); [] uses the default
%          all-items-outside-plus-first-m-in-list initial state.
%
% Returns the time grid TCACHE, the per-cache per-class hit/miss probability
% trajectories HITPROB_T/MISSPROB_T (ncaches x nclasses x nt), and the cache
% node indices CACHES (rows ordered as find(sn.nodetype==Cache)). XOCC carries
% the per-item, per-list occupancy trajectory of each cache and LAMCELL the
% isolated per-class per-item request rates that weight it, which is what a
% caller needs to resolve the hit probability by list rather than in total.
%
% HOLDMAP: optional holding-time MAP {D0,D1} of a random-environment stage.
%          When given, SOJOURN{c} is the average of cache c's transient
%          against that holding time, the drift's per-request time being
%          mapped to real time by the cache's total request rate LAM:
%          XBAR = int x(t) dF(t/LAM) / int dF(t/LAM) over OPTIONS.TIMESPAN,
%          with fields XBAR, WTOT (= the F mass over the window), HITPROB and
%          MISSPROB (1 x nclasses, NaN-free, zero for a class that does not
%          read the cache). The integral is carried by the ODE integrator
%          (CACHE_SOJOURN_ODE), so it does not depend on the output grid.
%          OPTIONS.TIMESPAN IS THEN REAL TIME, the stage horizon the holding
%          time is measured in: the drift is integrated over LAM*TIMESPAN in its
%          own per-request time and TCACHE is reported back in real time. Read
%          as request time, a window chosen to cover the holding time (a stage
%          horizon, or a few mean holding times) covered only 1/LAM of it, and
%          a busy cache's sojourn average was truncated to the early transient.
%
% The access graph of each cache (ch.accost) enters both the rate convergence
% and the transient drift, as it does in the steady SOLVER_FLD_CACHEQN_ANALYZER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
if nargin < 3
    x0cell = {};
end
if nargin < 4
    holdmap = [];
end
tspan = options.timespan;
if isinf(tspan(1))
    tspan(1) = 0;
end

% Converge the cache arrival rates and capture the per-cache isolated inputs.
[~, ~, ~, ~, ~, cacheinfo] = da_cacheqn(sn, @miss_isolated, @netsolve, options);
caches = cacheinfo.node;
ncaches = numel(caches);

tcache = [];
hitprob_t = [];
missprob_t = [];
arate = zeros(ncaches, K);
xocc = cell(1, ncaches);
lamcell = cell(1, ncaches);
sojourn = cell(1, ncaches);
for cIdx = 1:ncaches
    gamma = cacheinfo.gamma{cIdx};
    m = cacheinfo.m{cIdx};
    lam = cacheinfo.lambda_cache{cIdx};
    strat = cacheinfo.strat{cIdx};
    lamcell{cIdx} = lam;
    if cIdx <= numel(x0cell)
        x0 = x0cell{cIdx};
    else
        x0 = [];
    end

    % The holding-time clock in the drift's per-request time, when asked for.
    % The window is real time then, so the drift runs over LAM times it.
    wgen = [];
    tdrift = tspan;
    Lam = sum(sum(lam(:, :, 1)));
    if ~isempty(holdmap) && Lam > 0
        D0h = holdmap{1};
        tdrift = Lam * tspan;
        wgen = struct('A', D0h / Lam, 'c', -D0h * ones(size(D0h, 1), 1) / Lam);
        wgen.phi0 = map_pie(holdmap) * expm(D0h * tspan(1));
    end
    Rcost = cacheinfo.Rcost{cIdx};

    % Isolated-cache transient occupancy via the drift-based mean field. RANDOM(m)
    % and FIFO(m) share the steady state (Gast15 Thm 1) but NOT the transient, so
    % FIFO uses its own position-resolved drift; strict FIFO(m) likewise.
    % LRU/HLRU/CLIMB/QLRU have no drift-based transient.
    if strat == ReplacementStrategy.RR
        [~, ~, ~, ~, tc, ~, MU_t, xtraj, ~, ~, sj] = cache_miss_rmf(gamma, m, lam, tdrift, x0, Rcost, wgen);
        xocc{cIdx} = xtraj;
    elseif strat == ReplacementStrategy.FIFO
        [~, ~, ~, ~, tc, ~, MU_t, xtraj, sj] = cache_miss_fifo_rmf(gamma, m, lam, tdrift, x0, Rcost, wgen);
        xocc{cIdx} = xtraj;
    elseif strat == ReplacementStrategy.SFIFO
        [~, ~, ~, ~, tc, ~, MU_t, xtraj, sj] = cache_miss_sfifo_rmf(gamma, m, lam, tdrift, x0, Rcost, wgen);
        xocc{cIdx} = xtraj;
    else
        line_error(mfilename, sprintf(['Transient cache analysis is only ' ...
            'available for RANDOM(m)/FIFO(m) and strict FIFO(m) replacement via ' ...
            'a drift-based mean field; cache %d uses a strategy without one.'], cIdx));
    end

    if ~isequal(tdrift, tspan)
        tc = tc / Lam; % back to real time
    end
    if isempty(tcache)
        nt = numel(tc);
        tcache = tc(:)';
        hitprob_t = zeros(ncaches, K, nt);
        missprob_t = zeros(ncaches, K, nt);
    end

    % Per-user (per-class) arrival rate = sum over items of its isolated rate.
    u = size(lam, 1);
    for v = 1:u
        rowrate = sum(lam(v, :, 1));
        arate(cIdx, v) = rowrate;
        if rowrate > 0
            mp = MU_t(v, :) ./ rowrate;
            mp = max(0, min(1, mp));
            missprob_t(cIdx, v, :) = reshape(mp, 1, 1, []);
            hitprob_t(cIdx, v, :) = reshape(1 - mp, 1, 1, []);
        end
    end
    if ~isempty(sj)
        sojourn{cIdx} = struct('xbar', sj.xbar, 'wtot', sj.wtot, ...
            'hitprob', zeros(1, K), 'missprob', zeros(1, K));
        if ~isempty(sj.xbar)
            for v = 1:u
                if arate(cIdx, v) > 0
                    mpb = max(0, min(1, sj.MU(v) / arate(cIdx, v)));
                    sojourn{cIdx}.missprob(v) = mpb;
                    sojourn{cIdx}.hitprob(v) = 1 - mpb;
                end
            end
        end
    end
end

    function missrate = miss_isolated(gamma, m, lambda_cache, ch)
        % The rates the transient starts from are the steady analyzer's, access
        % graph included: FIFO(m) == RANDOM(m) only on the linear chain.
        if ch.replacestrat == ReplacementStrategy.RR
            [~, missrate] = cache_miss_rmf(gamma, m, lambda_cache, [], [], ch.accost);
        elseif ch.replacestrat == ReplacementStrategy.FIFO
            if isempty(cache_build_item_graphs(ch.accost, lambda_cache, size(lambda_cache, 2), length(m)))
                [~, missrate] = cache_miss_rmf(gamma, m, lambda_cache);
            else
                [~, missrate] = cache_miss_fifo_rmf(gamma, m, lambda_cache, [], [], ch.accost);
            end
        elseif ch.replacestrat == ReplacementStrategy.SFIFO
            [~, missrate] = cache_miss_sfifo_rmf(gamma, m, lambda_cache, [], [], ch.accost);
        else
            % LRU/HLRU/CLIMB/QLRU have no drift-based fluid model.
            line_error(mfilename, sprintf(['SolverFLD supports only ' ...
                'RANDOM(m)/FIFO(m) and strict FIFO(m) cache replacement; ' ...
                'strategy %d has no drift-based fluid model.'], ...
                double(ch.replacestrat)));
        end
    end

    function res = netsolve(snit)
        res = struct();
        fluid_options = options;
        fluid_options.method = 'matrix';
        fluid_options.init_sol = solver_fluid_initsol(snit, fluid_options);
        [res.QN, res.UN, res.RN, res.TN, res.xvec_iter, res.QNt, res.UNt, res.TNt, ~, res.t] = solver_fluid_matrix(snit, fluid_options);
        res.XN = zeros(1, K);
        for k = 1:K
            if snit.refstat(k) > 0
                res.XN(k) = res.TN(snit.refstat(k), k);
            end
        end
    end
end
