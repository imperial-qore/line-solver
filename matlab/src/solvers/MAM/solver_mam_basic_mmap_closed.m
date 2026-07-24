function [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap_closed(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_MAM_BASIC_MMAP_CLOSED(SN, OPTIONS)
%
% Closed-network wrapper around solver_mam_basic_mmap_inner. Drives a
% per-class bisection on the surrogate arrival rate LAMBDA so that the
% inner solver's queue lengths match the closed population SN.NJOBS.
% Mirrors the outer-loop structure of solver_mna_closed.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
M = sn.nstations;
S = 1./sn.rates;

% Per-class bisection bounds: upper = slowest non-INF station rate for that class
nonInfStations = find(sn.nservers < Inf);
lambda_lb = zeros(1,K);
lambda_ub = zeros(1,K);
for k=1:K
    rates_k = [];
    if ~isempty(nonInfStations)
        rates_k = sn.rates(nonInfStations, k);
        rates_k = rates_k(isfinite(rates_k) & rates_k > 0);
    end
    if isempty(rates_k)
        infStations = find(sn.nservers == Inf);
        rates_inf = sn.rates(infStations, k);
        rates_inf = rates_inf(isfinite(rates_inf) & rates_inf > 0);
        if isempty(rates_inf)
            lambda_ub(k) = 1;
        else
            lambda_ub(k) = max(rates_inf);
        end
    else
        lambda_ub(k) = min(rates_k);
    end
end

QNc = sn.njobs;
QNc(~isfinite(QNc)) = 0;  % open classes contribute 0; only closed populations gate convergence
QN_chain = zeros(1,K);

it_out = 0;
lambda = lambda_ub;
% Self-looping classes are pinned by the SLC clamp below; they must not
% contribute a (saturating) surrogate arrival stream to the inner algorithm.
lambda(sn.isslc) = 0;

inner_options = options;
inner_options.iter_max = max(20, ceil(options.iter_max/10));
inner_options.verbose = false;
% Cap MMAP phase truncation at 16 (lossless, avoids O(dim^3) waste);
% see _kb/06-solver-catalog.md for rationale
if ~isfield(inner_options.config, 'space_max') || inner_options.config.space_max > 16
    inner_options.config.space_max = 16;
end

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

% Last successful inner-algorithm outputs (for fallback if final trial diverges)
QN_last = QN; UN_last = UN; RN_last = RN;
TN_last = TN; CN_last = CN; XN_last = XN;
have_good = false;

bisect_tol = max(options.iter_tol, 1e-3);

while max(abs(QN_chain - QNc)) > bisect_tol && it_out < options.iter_max
    it_out = it_out + 1;
    if it_out > 1
        bracket_collapsed = true;
        for k=1:K
            if ~isfinite(QNc(k)) || QNc(k) == 0 || sn.isslc(k)
                continue;
            end
            if QN_chain(k) < QNc(k)
                lambda_lb(k) = lambda(k);
            else
                lambda_ub(k) = lambda(k);
            end
            lambda(k) = 0.5 * (lambda_lb(k) + lambda_ub(k));
            % Bisection can still refine class k only while its bracket is
            % wider than the precision floor below which LAMBDA cannot move
            % any reported metric.
            if (lambda_ub(k) - lambda_lb(k)) > GlobalConstants.FineTol * max(1, abs(lambda_ub(k)))
                bracket_collapsed = false;
            end
        end
        % Bracket-width stagnation break; see _kb/06-solver-catalog.md for rationale
        if bracket_collapsed
            it_out = it_out - 1;
            break;
        end
    end

    try
        [QN, UN, RN, TN, CN, XN, ~] = solver_mam_basic_mmap_inner(sn, inner_options, lambda);
        algorithm_ok = true;
    catch
        % Inner algorithm diverged (typically MMAPPH1FCFS / lyap NaN under
        % saturation). Treat all chains as overloaded so the bisection
        % drops lambda on its next step.
        algorithm_ok = false;
    end

    if algorithm_ok
        % SLC clamp: all jobs at refstat for self-looping classes
        for k=1:K
            if sn.isslc(k)
                QN(:,k) = 0;
                QN(sn.refstat(k), k) = sn.njobs(k);
            end
        end
        QN_chain = sum(QN, 1);
        QN_chain(isnan(QN_chain) | isinf(QN_chain)) = 1/GlobalConstants.FineTol;
        QN_last = QN; UN_last = UN; RN_last = RN;
        TN_last = TN; CN_last = CN; XN_last = XN;
        have_good = true;
    else
        QN_chain = ones(1,K) * (1/GlobalConstants.FineTol);
    end
end

% If the last trial diverged, fall back to the most recent successful one
if ~algorithm_ok && have_good
    QN = QN_last; UN = UN_last; RN = RN_last;
    TN = TN_last; CN = CN_last; XN = XN_last;
end

% Final SLC pass: pin throughput/utilisation at refstat (mirrors solver_mna_closed)
for k=1:K
    if sn.isslc(k)
        QN(:,k) = 0;
        ist = sn.refstat(k);
        QN(ist, k) = sn.njobs(k);
        TN(ist, k) = sn.njobs(k) * sn.rates(ist, k);
        if TN(ist, k) > 0
            RN(ist, k) = QN(ist, k) / TN(ist, k);
        else
            RN(ist, k) = 0;
        end
        UN(ist, k) = S(ist, k) * TN(ist, k);
    end
end

% Population redistribution within chain (matches solver_mna_closed:323-328)
for c=1:sn.nchains
    inchain = sn.inchain{c};
    if isfinite(sn.njobs(c))
        sumQ = sum(sum(QN(:,inchain)));
        if sumQ > 0
            QN(:,inchain) = sn.njobs(c) .* QN(:,inchain) / sumQ;
        end
    end
end

% Delay/INF utilisation = mean number of jobs (matches solver_mna_closed)
for ist=1:sn.nstations
    if sn.sched(ist) == SchedStrategy.INF
        UN(ist,:) = QN(ist,:);
    end
end

CN = sum(RN, 1);
QN(isnan(QN)) = 0;
UN(isnan(UN)) = 0;
RN(isnan(RN)) = 0;
TN(isnan(TN)) = 0;
CN(isnan(CN)) = 0;
totiter = it_out;
end
