function result = sim_firquest(Y, p, alpha, options)
% OA_FIRQUEST Fixed-sample-size quantile interval from independent replications.
%
% RESULT = OA_FIRQUEST(Y, P) returns a point estimate and a 95% confidence
% interval for the P-quantile of the steady-state marginal distribution, from R
% independent replications of equal length. Y is an N-by-R matrix whose columns
% are the replicate sample paths, or a cell array of R equal-length vectors.
%
% RESULT = OA_FIRQUEST(Y, P, ALPHA) uses nominal coverage 1-ALPHA.
%
% RESULT = OA_FIRQUEST(Y, P, ALPHA, OPTIONS) overrides the procedure constants;
% see OA_FQUEST for the recognized fields. Two defaults differ from FQUEST:
% b0 = 25 rather than 50, and s is chosen from R rather than fixed, because the
% stage tests act on the R*B pooled statistics and so need fewer batches per
% replication. Supplying s explicitly overrides that choice.
%
% The procedure is FIRQUEST, the replicated counterpart of FQUEST. It differs
% from OA_FQUEST in four places:
%
%   The warmup randomness test runs independently on each replicate path, and
%   the batch size it settles on may differ between replications.
%
%   Truncation removes the largest of those batch sizes from the front of every
%   replication, not just from one path. This is more aggressive than FQUEST on
%   purpose: an untruncated transient common to all replications biases every
%   replicate estimate the same way, and averaging cannot remove it.
%
%   The four stage tests act on the R*B signed areas and R*B replicate batched
%   quantile estimators pooled across replications, with the same B in every
%   replication and at least one batch from each.
%
%   The delivered interval is
%     ytilde_p(N*) +- t_{1-alpha/2, 2Rb-1} sqrt(Vtilde_p(w;R,b,m)/N*),
%   N* = R*B*M, with the pooled combined variance-parameter estimator
%     A_p(w;R,b,m)     = (Rb)^{-1} sum_j A_p(w;j,m)^2
%     Ntilde_p(R,b,m)  = m (Rb-1)^{-1} sum_j (yhat_p(j,m) - ytilde_p(N*))^2
%     Vtilde_p         = [Rb A_p + (Rb-1) Ntilde_p] / (2Rb-1).
%   The heuristic fallback drops FQUEST's residual-autocorrelation correction,
%   since the pooled batch quantiles come from independent paths.
%
% The default batch counts s as a function of R are the article's:
%   R = 2 -> [14 11 8 5], R = 3 -> [10 8 6 4], R = 4 -> [6 5 4 3],
%   5..9 -> [5 4 3 2], 10..16 -> [4 3 2 1], 17..22 -> [3 2 1],
%   23..32 -> [2 1], R >= 33 -> [1].
%
% Returns a struct with the fields OA_FQUEST returns, plus
%   R           - Number of replications
%   truncated   - Observations deleted from the front of every replication
% and with b and m the per-replication batch count and batch size, so that
% n = R*b*m.
%
% Independent replications shorten the correlation the estimator has to fight,
% and they parallelize, but they reintroduce initialization bias in every path,
% so a short run length per replication is worse here than in OA_FQUEST. The
% article reports slight undercoverage at P = 0.99 when the total sample is
% under 500000, down to 90.8%.
%
% Examples:
%   Y = exprnd(1, 40000, 5);
%   r = sim_firquest(Y, 0.9, 0.05);
%   [r.lower r.estimate r.upper]     % brackets -log(0.1) = 2.3026
%
% Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec, A. C. Mokashi,
% J. R. Wilson, "A Fixed-Sample-Size Procedure for Estimating Steady-State
% Quantiles Based on Independent Replications", Proc. Winter Simulation
% Conference, 2025.
%
% See also OA_FQUEST, OA_STS_QUANTILE_AREAS, OA_VONNEUMANN, OA_SHAPIROWILK
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 4
    options = struct();
end

if iscell(Y)
    lens = cellfun(@numel, Y);
    if isempty(Y) || any(lens ~= lens(1))
        line_error(mfilename, ...
            'the replicate paths must be nonempty and of equal length');
    end
    Ym = zeros(lens(1), numel(Y));
    for r = 1:numel(Y)
        col = Y{r};
        Ym(:, r) = col(:);
    end
else
    Ym = Y;
    if isvector(Ym)
        line_error(mfilename, ...
            'at least 2 replications are required, use sim_fquest for a single path');
    end
end

[n, R] = size(Ym);
if R < 2
    line_error(mfilename, ...
        'at least 2 replications are required, use sim_fquest for a single path');
end
if any(~isfinite(Ym(:)))
    line_error(mfilename, 'the sample paths must be finite');
end

if ~isfield(options, 'b0')
    options.b0 = 25;
end
if ~isfield(options, 's')
    options.s = sim_firquest_batchcounts(R);
end
opt = sim_quest_options(options);

if ~isscalar(p) || ~isreal(p) || p <= 0 || p >= 1
    line_error(mfilename, 'p must be a real scalar in (0,1)');
end
if ~isscalar(alpha) || ~isreal(alpha) || alpha <= 0 || alpha >= 1
    line_error(mfilename, 'alpha must be a real scalar in (0,1)');
end
if R * opt.s(end) < 3
    line_error(mfilename, ...
        'R*min(s) = %d pooled batches is below the 3 the stage tests need', ...
        R * opt.s(end));
end

warnings = {};

% ---- warmup: one randomness loop per replicate path
b = opt.b0;
mStart = opt.m0;
if n < b * mStart
    mStart = floor(n / b);
end
if mStart < 1
    line_error(mfilename, ...
        'each replication holds %d observations, too few for b0 = %d batches', n, b);
end

mMax = 0;
failed = false;
for r = 1:R
    m = mStart;
    ell = 1;
    atMax = false;
    passed = false;
    while true
        stats = sim_sts_quantile_areas(Ym(1:b * m, r), b, m, p, opt.weight);
        sig = opt.beta * exp(-opt.eta * (ell - 1)^opt.theta);
        if ~sim_vonneumann(stats.areas, sig).reject
            passed = true;
            break
        end
        if atMax
            break
        end
        ell = ell + 1;
        mNext = round(m * sqrt(2));
        if n < b * mNext && mNext ~= floor(n / b)
            m = floor(n / b);
        else
            m = mNext;
            if n < b * m
                m = floor(n / b);
                atMax = true;
            end
        end
        if m < 1
            break
        end
    end
    failed = failed || ~passed;
    mMax = max(mMax, m);
end
if failed
    warnings{end + 1} = ['the warmup randomness test could not be passed in ' ...
        'every replication, the replicate paths are too short'];
end

% ---- truncation: delete the longest warmup batch from every replication
truncated = mMax;
if truncated >= n
    line_error(mfilename, ...
        'the warmup batch size %d exhausts the replication length %d', truncated, n);
end
Yt = Ym(truncated + 1:end, :);
nstarRep = size(Yt, 1);

% ---- batch-count selection on the pooled statistics
v = 1;
b = opt.s(v);
m = floor(nstarRep / b);
ok = true;
pooled = [];
for stage = 1:4
    while true
        if m < 1
            ok = false;
            break
        end
        pooled = sim_firquest_pool(Yt, b, m, p, opt.weight);
        if stage <= 2
            sample = pooled.areas;
        else
            sample = pooled.bqe;
        end
        if mod(stage, 2) == 1
            reject = sim_vonneumann(sample, opt.beta).reject;
        else
            reject = sim_shapirowilk(sample, opt.beta).reject;
        end
        if ~reject
            break
        end
        v = v + 1;
        if v > numel(opt.s)
            ok = false;
            break
        end
        b = opt.s(v);
        m = floor(nstarRep / b);
    end
    if ~ok
        break
    end
end

if isempty(pooled) || m < 1
    line_error(mfilename, ...
        'the replicate paths are too short to form %d batches each', opt.s(end));
end

nstar = pooled.n;
estimate = pooled.quantile;

if ok
    half = sim_tinv(1 - alpha / 2, 2 * R * b - 1) * sqrt(pooled.Vp / nstar);
    lower = estimate - half;
    upper = estimate + half;
    heuristic = false;
else
    warnings{end + 1} = sprintf(['a randomness or normality test failed at ' ...
        'b = %d per replication, the delivered interval is heuristic'], opt.s(end));
    heuristic = true;
    if opt.force
        [lower, upper] = sim_quest_heuristic_ci(pooled.bqe, estimate, pooled.Ap, ...
            pooled.Np, nstar, alpha, false);
        half = (upper - lower) / 2;
    else
        lower = NaN;
        upper = NaN;
        half = NaN;
    end
end

result = struct('estimate', estimate, 'lower', lower, 'upper', upper, ...
    'halfwidth', half, 'b', b, 'm', m, 'n', nstar, 'R', R, ...
    'truncated', truncated, 'Ap', pooled.Ap, 'Np', pooled.Np, 'Vp', pooled.Vp, ...
    'heuristic', heuristic, 'warnings', {warnings}, 'analyzer', 'sim_firquest');
end

function s = sim_firquest_batchcounts(R)
% Article default batch counts, chosen so that R*b pooled statistics remain
% enough to test while every replication still contributes at least one batch.
if R == 2
    s = [14 11 8 5];
elseif R == 3
    s = [10 8 6 4];
elseif R == 4
    s = [6 5 4 3];
elseif R < 10
    s = [5 4 3 2];
elseif R < 17
    s = [4 3 2 1];
elseif R < 23
    s = [3 2 1];
elseif R < 33
    s = [2 1];
else
    s = 1;
end
end

function pooled = sim_firquest_pool(Yt, b, m, p, weight)
% Pool the per-replication areas and batch quantiles, then form the replicated
% variance-parameter estimators over all R*b of them.
R = size(Yt, 2);
nstarRep = size(Yt, 1);
kept = Yt(nstarRep - b * m + 1:end, :);
areas = zeros(b, R);
bqe = zeros(b, R);
for r = 1:R
    stats = sim_sts_quantile_areas(kept(:, r), b, m, p, weight);
    areas(:, r) = stats.areas;
    bqe(:, r) = stats.bqe;
end

n = R * b * m;
allSorted = sort(kept(:));
quantile = allSorted(ceil(n * p));

K = R * b;
Ap = sum(areas(:).^2) / K;
Np = m * sum((bqe(:) - quantile).^2) / (K - 1);
Vp = (K * Ap + (K - 1) * Np) / (2 * K - 1);

pooled = struct('areas', areas(:), 'bqe', bqe(:), 'quantile', quantile, ...
    'Ap', Ap, 'Np', Np, 'Vp', Vp, 'b', b, 'm', m, 'n', n, 'R', R);
end
