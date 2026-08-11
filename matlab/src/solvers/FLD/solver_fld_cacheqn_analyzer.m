function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec_iter, hitprob, missprob, runtime, it] = solver_fld_cacheqn_analyzer(sn, options)
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

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
K = sn.nclasses;

[res, hitprob, missprob, it, sn] = da_cacheqn(sn, @miss_isolated, @netsolve, options);
QN = res.QN; UN = res.UN; RN = res.RN; TN = res.TN; XN = res.XN;
t = res.t; QNt = res.QNt; UNt = res.UNt; TNt = res.TNt; xvec_iter = res.xvec_iter;

% Compute CN
CN = zeros(1, K);
for k = 1:K
    if sn.refstat(k) > 0
        CN(k) = sn.njobs(k) ./ XN(k);
    end
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
        % Solve the queueing network using the fluid matrix method
        res = struct();
        fluid_options = options;
        fluid_options.method = 'matrix';
        fluid_options.init_sol = solver_fluid_initsol(snit, fluid_options);
        [res.QN, res.UN, res.RN, res.TN, res.xvec_iter, res.QNt, res.UNt, res.TNt, ~, res.t] = solver_fluid_matrix(snit, fluid_options);

        % Compute system throughputs
        res.XN = zeros(1, K);
        for k = 1:K
            if snit.refstat(k) > 0
                res.XN(k) = res.TN(snit.refstat(k), k);
            end
        end
    end
end
