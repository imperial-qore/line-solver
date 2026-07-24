function [tcache, hitprob_t, missprob_t, caches, arate, xocc] = solver_fld_cacheqn_tran(sn, options, x0cell)
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
% node indices CACHES (rows ordered as find(sn.nodetype==Cache)).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
if nargin < 3
    x0cell = {};
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
for cIdx = 1:ncaches
    gamma = cacheinfo.gamma{cIdx};
    m = cacheinfo.m{cIdx};
    lam = cacheinfo.lambda_cache{cIdx};
    strat = cacheinfo.strat{cIdx};
    if cIdx <= numel(x0cell)
        x0 = x0cell{cIdx};
    else
        x0 = [];
    end

    % Isolated-cache transient occupancy via the drift-based mean field. RANDOM(m)
    % and FIFO(m) share the steady state (Gast15 Thm 1) but NOT the transient, so
    % FIFO uses its own position-resolved drift; strict FIFO(m) likewise.
    % LRU/HLRU/CLIMB/QLRU have no drift-based transient.
    if strat == ReplacementStrategy.RR
        [~, ~, ~, ~, tc, ~, MU_t, xtraj] = cache_miss_rmf(gamma, m, lam, tspan, x0);
        xocc{cIdx} = xtraj;
    elseif strat == ReplacementStrategy.FIFO
        [~, ~, ~, ~, tc, ~, MU_t, xtraj] = cache_miss_fifo_rmf(gamma, m, lam, tspan, x0);
        xocc{cIdx} = xtraj;
    elseif strat == ReplacementStrategy.SFIFO
        [~, ~, ~, ~, tc, ~, MU_t, xtraj] = cache_miss_sfifo_rmf(gamma, m, lam, tspan, x0);
        xocc{cIdx} = xtraj;
    else
        line_error(mfilename, sprintf(['Transient cache analysis is only ' ...
            'available for RANDOM(m)/FIFO(m) and strict FIFO(m) replacement via ' ...
            'a drift-based mean field; cache %d uses a strategy without one.'], cIdx));
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
end

    function missrate = miss_isolated(gamma, m, lambda_cache, ch)
        if ch.replacestrat == ReplacementStrategy.RR || ...
                ch.replacestrat == ReplacementStrategy.FIFO
            [~, missrate] = cache_miss_rmf(gamma, m, lambda_cache);
        elseif ch.replacestrat == ReplacementStrategy.SFIFO
            [~, missrate] = cache_miss_sfifo_rmf(gamma, m, lambda_cache);
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
