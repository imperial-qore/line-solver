function result = sim_fquest(Y, p, alpha, options)
% OA_FQUEST Fixed-sample-size confidence interval for a steady-state quantile.
%
% RESULT = OA_FQUEST(Y, P) returns a point estimate and a 95% confidence
% interval for the P-quantile of the steady-state marginal distribution of the
% simulation output process whose single sample path is Y. Y has arbitrary
% fixed length; no sequential control of the run length is needed.
%
% RESULT = OA_FQUEST(Y, P, ALPHA) uses nominal coverage 1-ALPHA.
%
% RESULT = OA_FQUEST(Y, P, ALPHA, OPTIONS) overrides the procedure constants.
% OPTIONS is a struct whose recognized fields are
%   b0       - initial batch count for the warmup stage (default 50)
%   m0       - initial batch size for the warmup stage (default 500)
%   s        - descending batch counts for the test stages (default [32 24 16 10])
%   beta     - significance level of the stage tests (default 0.30)
%   eta      - decay coefficient of the warmup significance (default 0.2)
%   theta    - decay exponent of the warmup significance (default 2.3)
%   weight   - constant STS weight function (default sqrt(12))
%   force    - deliver a heuristic interval when a test fails (default true)
%
% The procedure is FQUEST. It has four blocks:
%
%   Warmup. Starting from B = b0 and M = m0 it computes the B signed STS areas
%   of the batched quantile process and tests them for randomness with von
%   Neumann's ratio at the decaying significance beta*exp(-eta*(l-1)^theta) on
%   iteration l, growing M by sqrt(2) whenever the test rejects. Passing the
%   test means the areas are approximately independent, so any initialization
%   bias is confined to the first batch.
%
%   Truncation. The first batch is deleted, which is the entire warmup
%   treatment; there is no separate transient detector.
%
%   Batch-count selection. With B stepping down through s and M = floor(N*/B),
%   four tests must pass in order: von Neumann and Shapiro-Wilk on the signed
%   areas, then von Neumann and Shapiro-Wilk on the batched quantile
%   estimators. These check the asymptotic properties the interval rests on,
%   namely that both the areas and the batch quantiles behave like independent
%   normal variates. B only ever decreases, and a failure at B = 10 ends the
%   stage.
%
%   Delivery. When all four tests pass the interval is
%     ytilde_p(n*) +- t_{1-alpha/2, 2b-1} sqrt(V_p(w;b,m)/n*)
%   with V_p the combined variance-parameter estimator of
%   OA_STS_QUANTILE_AREAS. Otherwise the sample was too small, RESULT.heuristic
%   is true, and with OPTIONS.force the interval returned is the union of the
%   wider of the two single-component intervals and Willink's skewness- and
%   correlation-adjusted asymmetric interval, both built from the batched
%   quantile estimators.
%
% Returns a struct with fields:
%   estimate    - Full-sample empirical P-quantile of the truncated path
%   lower,upper - Confidence interval endpoints
%   halfwidth   - Half-width, (upper-lower)/2, which the asymmetric heuristic
%                 interval attains only on average
%   b, m, n     - Final batch count, batch size and number of observations used
%   truncated   - Number of observations deleted from the front of Y
%   Ap, Np, Vp  - The three variance-parameter estimators at the final b and m
%   heuristic   - true when a test failed and the interval is not asymptotically
%                 justified
%   warnings    - Cell array of diagnostic strings, empty on a clean run
%   analyzer    - Identifier string
%
% Coverage was measured on the article's own test bed, the waiting-time process
% of an M/M/1 queue with lambda = 0.8, mu = 1 started with 113 jobs in system,
% over 500 independent replications at N = 200000, giving a standard error near
% 1%: 95.2% at p = 0.5, 96.2% at p = 0.9 and 95.6% at p = 0.99 against a nominal
% 95%. The delivered half-width exceeds the empirically needed one by factors of
% 1.16, 1.24 and 1.99 respectively, so the interval is conservative and
% increasingly so into the tail, consistent with the half-widths the article
% reports. On i.i.d. Exp(1) data, where sigma_p^2 = p(1-p)/f(y_p)^2 is exact,
% both A_p and N_p are unbiased to within 7%.
%
% Two properties are worth knowing. Coverage is not monotone in the sample size
% over this range: a larger sample passes the stage tests more often and so
% reaches the conservative fallback less. And a substantial fraction of runs
% takes that fallback at all, 22% to 65% here and rising with p, so a delivered
% interval may well be the heuristic one; RESULT.heuristic says which. Fewer
% than about 100 replications cannot resolve a two-point difference in coverage,
% so do not read a small experiment as a defect.
%
% Applicability is a condition on the output process, not on the model that
% produced it. The theory needs geometric moment contraction (Wu 2005), which
% holds for ARMA series, a broad class of short-range-dependent linear and
% nonlinear processes, many Markov chains, and was proved for M/M/1 and
% non-heavy-tailed G/G/1 waiting times by Dingec et al. (2022); a density that is
% positive and differentiable at the quantile of interest; short-range dependence
% and an FCLT for the indicator process. M/M/1 is only the validation bed, chosen
% because its exact quantiles are known.
%
% Two practical exclusions follow. Do NOT use this on integer-valued output such
% as a queue length: the marginal has no density, the density-regularity
% condition fails, and the batched quantile has no Bahadur representation. Use it
% on continuous output, that is response, waiting and sojourn times. And
% heavy-tailed service, which can break geometric moment contraction and induce
% long-range dependence, is outside the theory.
%
% Examples:
%   Y = exprnd(1, 200000, 1);
%   r = sim_fquest(Y, 0.9, 0.05);
%   [r.lower r.estimate r.upper]     % brackets -log(0.1) = 2.3026
%
% Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec, A. C. Mokashi,
% J. R. Wilson, "A Fixed-Sample-Size Method for Estimating Steady-State
% Quantiles", Proc. Winter Simulation Conference, 2023.
%
% See also OA_FIRQUEST, OA_STS_QUANTILE_AREAS, OA_VONNEUMANN, OA_SHAPIROWILK
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 4
    options = struct();
end
opt = sim_quest_options(options);

if ~isscalar(p) || ~isreal(p) || p <= 0 || p >= 1
    line_error(mfilename, 'p must be a real scalar in (0,1)');
end
if ~isscalar(alpha) || ~isreal(alpha) || alpha <= 0 || alpha >= 1
    line_error(mfilename, 'alpha must be a real scalar in (0,1)');
end

if opt.s(end) < 3
    line_error(mfilename, ...
        'min(s) = %d, but the stage tests need at least 3 batches', opt.s(end));
end

Y = Y(:);
N = numel(Y);
if any(~isfinite(Y))
    line_error(mfilename, 'the sample path must be finite');
end
if N < opt.s(end) * 2
    line_error(mfilename, ...
        'the sample path holds %d observations, at least %d are needed', ...
        N, opt.s(end) * 2);
end

warnings = {};

% ---- warmup: grow the batch size until the signed areas look random
b = opt.b0;
m = opt.m0;
if N < b * m
    m = floor(N / b);
end
if m < 1
    line_error(mfilename, ...
        'the sample path is too short for the initial batch count b0 = %d', b);
end

ell = 1;
atMax = false;
passed = false;
while true
    stats = sim_sts_quantile_areas(Y(1:b * m), b, m, p, opt.weight);
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
    if N < b * mNext && mNext ~= floor(N / b)
        m = floor(N / b);
    else
        m = mNext;
        if N < b * m
            m = floor(N / b);
            atMax = true;
        end
    end
    if m < 1
        break
    end
end
if ~passed
    warnings{end + 1} = ['the warmup randomness test could not be passed at ' ...
        'the largest admissible batch size, the sample path is too short'];
end

% ---- truncation: delete the first batch
truncated = m;
Yt = Y(truncated + 1:end);
Nstar = numel(Yt);

% ---- batch-count selection: four tests in order, b only decreases
v = 1;
b = opt.s(v);
m = floor(Nstar / b);
ok = true;
stats = [];
for stage = 1:4
    while true
        if m < 1
            ok = false;
            break
        end
        stats = sim_sts_quantile_areas(Yt(Nstar - b * m + 1:end), b, m, p, opt.weight);
        if stage <= 2
            sample = stats.areas;
        else
            sample = stats.bqe;
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
        m = floor(Nstar / b);
    end
    if ~ok
        break
    end
end

if isempty(stats) || m < 1
    line_error(mfilename, ...
        'the sample path is too short to form %d batches', opt.s(end));
end

nstar = stats.n;
estimate = stats.quantile;

if ok
    half = sim_tinv(1 - alpha / 2, 2 * b - 1) * sqrt(stats.Vp / nstar);
    lower = estimate - half;
    upper = estimate + half;
    heuristic = false;
else
    warnings{end + 1} = sprintf(['a randomness or normality test failed at ' ...
        'b = %d, the delivered interval is heuristic'], opt.s(end));
    heuristic = true;
    if opt.force
        [lower, upper] = sim_quest_heuristic_ci(stats.bqe, estimate, stats.Ap, ...
            stats.Np, nstar, alpha, true);
        half = (upper - lower) / 2;
    else
        lower = NaN;
        upper = NaN;
        half = NaN;
    end
end

result = struct('estimate', estimate, 'lower', lower, 'upper', upper, ...
    'halfwidth', half, 'b', b, 'm', m, 'n', nstar, 'truncated', truncated, ...
    'Ap', stats.Ap, 'Np', stats.Np, 'Vp', stats.Vp, 'heuristic', heuristic, ...
    'warnings', {warnings}, 'analyzer', 'sim_fquest');
end
