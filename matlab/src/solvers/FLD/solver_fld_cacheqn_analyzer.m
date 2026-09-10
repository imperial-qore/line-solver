function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec_iter, hitprob, missprob, runtime, it, momentResults] = solver_fld_cacheqn_analyzer(sn, options)
% SOLVER_FLD_CACHEQN_ANALYZER Fluid solver for integrated caching-queueing networks
%
% Delegates the decomposition-aggregation alternation between the isolated
% caches and the fluid ODE solution of the surrounding queueing network to
% da_cacheqn, until the cache arrival rates converge. RANDOM(m) (RR) and
% FIFO(m) use the refined mean field (cache_miss_rmf, 1/N-accurate; Gast15
% Thm 1: pi_FIFO(m)=pi_RAND(m)); strict FIFO(m) uses its own position-resolved
% mean field (cache_miss_sfifo_rmf). LRU/HLRU/CLIMB/QLRU have no drift-based
% fluid model and are rejected (the FPI characteristic-time approximation is
% not a fluid method). Use SolverNC/SolverMVA or SolverLDES for those.
%
% OPTIONS.METHOD SELECTS THE QUEUEING LAYER, not the cache one. 'rmf' solves
% the surrounding network with the first-order matrix method, which is the
% historical behaviour; 'minnormal' solves it with the second-order moment
% closure instead, so a cache model reaches the same E[min(X,c)] treatment as
% any other model and returns a covariance. The cache layer is the refined
% mean field either way -- it has no first-order alternative here -- and under
% 'minnormal' its own linear noise covariance is reported alongside, which is
% the second moment of the item occupancy.
%
% MOMENTRESULTS is empty unless the moment closure ran. It then carries the
% queueing fields of SOLVER_FLUID_MOMENTS (Sigma, QVar, QStd, sigma2,
% stationBlock, classBlock) plus CACHE, a struct array with one entry per cache
% node: NODE, PI0 (per-item miss probability), SIGMA (occupancy covariance,
% item-major), PI0VAR (per-item miss-indicator variance) and MISSPROBVAR (per
% class, the delta-method variance of the miss probability that class sees).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
K = sn.nclasses;
useMoments = any(strcmp(options.method, {'minnormal','fluid.minnormal'}));
lastMoments = [];

[res, hitprob, missprob, it, sn, cacheinfo] = da_cacheqn(sn, @miss_isolated, @netsolve, options);
QN = res.QN; UN = res.UN; RN = res.RN; TN = res.TN; XN = res.XN;
t = res.t; QNt = res.QNt; UNt = res.UNt; TNt = res.TNt; xvec_iter = res.xvec_iter;

% Compute CN
CN = zeros(1, K);
for k = 1:K
    if sn.refstat(k) > 0
        CN(k) = sn.njobs(k) ./ XN(k);
    end
end

% The moment report: the queueing fields as any other model returns them, plus
% the cache occupancy covariance evaluated at the converged isolated-cache
% inputs. The cache is re-solved once here rather than inside the sweep because
% only the FINAL arrival rates define the fixed point the covariance linearises
% about.
momentResults = [];
if useMoments && ~isempty(lastMoments)
    momentResults = lastMoments;
    momentResults.cache = local_cache_moments(cacheinfo, K);
end

% xvec_iter is already a cell array from solver_fluid_matrix
runtime = toc(T0);

    function missrate = miss_isolated(gamma, m, lambda_cache, ch)
        % RR/FIFO/SFIFO honour a custom access graph (accost) via their general
        % drift; the linear default keeps the refined RAND (FIFO) / linear
        % position-resolved (SFIFO) path. FIFO(m) == RANDOM(m) only for the
        % linear chain, so a non-linear graph uses its own position-resolved drift.
        nonlinear = ~accost_is_linear(ch.accost, length(m));
        if ch.replacestrat == ReplacementStrategy.RR
            [~, missrate] = cache_miss_rmf(gamma, m, lambda_cache, [], [], ch.accost);
        elseif ch.replacestrat == ReplacementStrategy.FIFO
            if nonlinear
                [~, missrate] = cache_miss_fifo_rmf(gamma, m, lambda_cache, [], [], ch.accost);
            else
                [~, missrate] = cache_miss_rmf(gamma, m, lambda_cache);
            end
        elseif ch.replacestrat == ReplacementStrategy.SFIFO
            [~, missrate] = cache_miss_sfifo_rmf(gamma, m, lambda_cache, [], [], ch.accost);
        else
            % LRU/HLRU/CLIMB/QLRU have no drift-based fluid model; the FPI
            % (characteristic-time) approximation is not a fluid method. Refuse
            % rather than substitute a non-fluid fixed point.
            line_error(mfilename, sprintf(['SolverFLD supports only ' ...
                'RANDOM(m)/FIFO(m) (refined mean field) and strict FIFO(m) ' ...
                '(position-resolved mean field) cache replacement; strategy ' ...
                '%d has no drift-based fluid model. Use SolverNC/SolverMVA or ' ...
                'SolverLDES for this cache.'], double(ch.replacestrat)));
        end
    end

    function cachemom = local_cache_moments(ci, nclasses)
        % Second moment of each cache, at the converged isolated-cache inputs.
        %
        % CACHE_MISS_RMF returns the stationary covariance of the item
        % occupancy under the linear noise approximation. The per-item miss
        % indicator is coordinate (i, list 0), so its variance is the leading
        % n_items block of the diagonal; the miss probability a class sees is
        % the popularity-weighted sum of those indicators, hence a linear
        % functional whose variance is w'*W00*w. Only RR/FIFO on the linear
        % access chain have the refined path that produces W, so a cache
        % without one reports empty rather than a fabricated zero.
        ncache = numel(ci.node);
        cachemom = struct('node', {}, 'pi0', {}, 'Sigma', {}, 'pi0Var', {}, 'missProbVar', {});
        for c = 1:ncache
            strat = ci.strat{c};
            if ~(strat == ReplacementStrategy.RR || strat == ReplacementStrategy.FIFO)
                continue
            end
            lam = ci.lambda_cache{c};
            if ~accost_is_linear(ci.Rcost{c}, length(ci.m{c}))
                continue
            end
            [~, ~, ~, pi0, ~, ~, ~, ~, ~, W] = cache_miss_rmf(ci.gamma{c}, ci.m{c}, lam);
            if isempty(W)
                continue
            end
            nitems = numel(pi0);
            W00 = W(1:nitems, 1:nitems);
            entry = struct();
            entry.node = ci.node(c);
            entry.pi0 = pi0(:);
            entry.Sigma = W;
            entry.pi0Var = max(0, diag(W00));
            entry.missProbVar = zeros(1, nclasses);
            for r = 1:min(nclasses, size(lam,1))
                w = lam(r,:,1);
                w(~isfinite(w)) = 0;
                tot = sum(w);
                if tot <= 0
                    continue
                end
                w = w(:) / tot;
                entry.missProbVar(r) = max(0, w' * W00 * w);
            end
            cachemom(end+1) = entry; %#ok<AGROW>
        end
    end

    function tf = accost_is_linear(accost, h)
        % True when every per-(user,item) access graph is the linear chain
        % (miss -> list 1, hit in list a -> list a+1, self-loop on top list).
        if isempty(accost)
            tf = true; return;
        end
        lin = zeros(h+1, h+1); lin(1,2) = 1;
        for a = 1:(h-1), lin(a+1, a+2) = 1; end
        lin(h+1, h+1) = 1;
        tf = true;
        [uu, nn] = size(accost);
        for v = 1:uu
            for kk = 1:nn
                g = accost{v, kk};
                if isempty(g), continue; end
                if ~all(all(abs(g - lin) < 1e-9))
                    tf = false; return;
                end
            end
        end
    end

    function res = netsolve(snit)
        % Solve the queueing network. The caches are already relabeled as class
        % switches by DA_CACHEQN, so SNIT is a plain queueing network and the
        % moment closure applies to it unchanged.
        res = struct();
        fluid_options = options;
        if useMoments
            fluid_options.method = 'minnormal';
            fluid_options.init_sol = solver_fluid_initsol(snit, fluid_options);
            [res.QN, res.UN, res.RN, res.TN, res.xvec_iter, res.QNt, res.UNt, res.TNt, ~, res.t, ~, ~, lastMoments] = ...
                solver_fluid_moments(snit, fluid_options);
        else
            fluid_options.method = 'matrix';
            fluid_options.init_sol = solver_fluid_initsol(snit, fluid_options);
            [res.QN, res.UN, res.RN, res.TN, res.xvec_iter, res.QNt, res.UNt, res.TNt, ~, res.t] = solver_fluid_matrix(snit, fluid_options);
        end

        % Compute system throughputs
        res.XN = zeros(1, K);
        for k = 1:K
            if snit.refstat(k) > 0
                res.XN(k) = res.TN(snit.refstat(k), k);
            end
        end
    end
end
