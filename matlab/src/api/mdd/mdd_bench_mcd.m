function T = mdd_bench_mcd(Ks, Ns, options)
% T = MDD_BENCH_MCD(KS, NS)
% T = MDD_BENCH_MCD(KS, NS, OPTIONS)
% Cost profile of the Miner-Ciardo-Donatelli aggregation (mdd_mcd) against the
% exact CTMC solution of the same model, over a grid of station counts KS and
% populations NS on single-class closed cyclic exponential networks.
%
% Three costs are separated because they are not interchangeable:
%   t_reach : building the exact reachable set into the MDD. Both methods need
%             it -- mdd_mcd reads the diagram, so it can never beat this term.
%   t_exact : explicit generator assembly + ctmc_solve, given the diagram.
%   t_mcd   : descriptor + the K coupled level-CTMC fixed point, given the
%             diagram.
% The crossover reported is therefore t_exact/t_mcd (the solve-phase speedup,
% the quantity the method is actually about) alongside the end-to-end ratio
% (t_reach+t_exact)/(t_reach+t_mcd), which is what a user observes.
%
% Storage is compared as |S| (exact stationary vector) against sum_k |M_k| (the
% K level vectors that mdd_mcd stores instead).
%
% -- Input
% KS      : vector of station counts (default [3 4 5 6])
% NS      : vector of populations   (default [4 8 12 16 20])
% OPTIONS : struct, fields
%             maxstates  - skip the exact solve above this |S| (default 3e4)
%             ctmcstates - also time LINE's SolverCTMC below this |S|
%                          (default 0 = never; SolverCTMC re-derives the state
%                          space itself and is the end-user baseline)
%             mcdopts    - options struct forwarded to mdd_mcd
%             reps       - repetitions per phase, median reported (default 3)
%             ctmcmethod - forced ctmc_solve method for the exact baseline;
%                          'direct' keeps the sparse factorization above
%                          GMRES_MIN_STATES (6000), where the default path
%                          switches to GMRES and its cost stops tracking |S|
%             verbose    - print the table as it is produced (default true)
% -- Output
% T       : table with one row per (K,N) case
%
% Temporary developer harness (not part of the public API).
%
% See also: mdd_mcd, mdd_closedqn, mdd_profile.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(Ks), Ks = [3 4 5 6]; end
if nargin < 2 || isempty(Ns), Ns = [4 8 12 16 20]; end
if nargin < 3 || isempty(options), options = struct(); end
if ~isfield(options, 'maxstates'),  options.maxstates = 3e4;  end
if ~isfield(options, 'ctmcstates'), options.ctmcstates = 0;   end
if ~isfield(options, 'mcdopts'),    options.mcdopts = struct(); end
if ~isfield(options, 'reps'),       options.reps = 3;         end
if ~isfield(options, 'ctmcmethod'), options.ctmcmethod = ''; end
if ~isfield(options, 'verbose'),    options.verbose = true;   end

rows = struct('K', {}, 'N', {}, 'nS', {}, 'sumMk', {}, 'maxMk', {}, 'nnzQ', {}, ...
    't_reach', {}, 't_exact', {}, 't_mcd', {}, 't_solverctmc', {}, ...
    'iters', {}, 'errQLen', {}, 'speedup', {}, 'speedup_e2e', {});

if options.verbose
    line_printf('\n%3s %4s %9s %9s %8s %9s %9s %9s %9s %6s %10s %9s %9s\n', ...
        'K', 'N', '|S|', 'sum|Mk|', 'max|Mk|', 'nnz(Q)', 't_reach', ...
        't_exact', 't_mcd', 'iters', 'errQLen', 'speedup', 'e2e');
end

for K = Ks(:)'
    for N = Ns(:)'
        nS = i_numstates(N, K);
        if nS > options.maxstates
            if options.verbose
                line_printf('%3d %4d %9d %9s %8s %9s %9s %9s %9s %6s %10s %9s %9s  (skipped: |S| > maxstates)\n', ...
                    K, N, nS, '-', '-', '-', '-', '-', '-', '-', '-', '-', '-');
            end
            continue
        end

        [mu, P, servers] = i_case(K);

        % ---- shared phase: exact reachable set into the diagram
        domain = (N + 1) * ones(1, K);
        init = zeros(1, K); init(1) = N;
        [ii, jj] = find(P .* (1 - eye(K)) > 0);
        treps = zeros(1, options.reps);
        for r = 1:options.reps
            t0 = tic;
            mdd = mdd_reachset(domain, init, @(s) i_next(s, ii, jj));
            treps(r) = toc(t0);
        end
        t_reach = median(treps);
        mdds = mdd.toStruct();

        % ---- exact solve, reusing the diagram so only the solve phase is timed
        eopt = struct('verbose', false, 'mdd', mdd);
        if ~isempty(options.ctmcmethod), eopt.ctmcmethod = options.ctmcmethod; end
        for r = 1:options.reps
            oute = mdd_closedqn(mu, P, servers, N, eopt);
            treps(r) = oute.times.gen + oute.times.solve + oute.times.metrics;
        end
        t_exact = median(treps);
        nnzQ = nnz(oute.Q);

        % ---- aggregation given the same diagram
        for r = 1:options.reps
            t0 = tic;
            desc = mdd_descriptor(mu, P, servers, N);
            outm = mdd_mcd(mdds, desc, options.mcdopts);
            treps(r) = toc(t0);
        end
        t_mcd = median(treps);

        errQ = max(abs(oute.QLen - outm.QLen));

        % ---- optional end-user baseline: LINE SolverCTMC on the same model
        t_ctmc = NaN;
        if nS <= options.ctmcstates
            t0 = tic;
            i_solverctmc(mu, servers, N);
            t_ctmc = toc(t0);
        end

        r = struct('K', K, 'N', N, 'nS', nS, ...
            'sumMk', sum(outm.levelSizes), 'maxMk', max(outm.levelSizes), ...
            'nnzQ', nnzQ, 't_reach', t_reach, 't_exact', t_exact, ...
            't_mcd', t_mcd, 't_solverctmc', t_ctmc, 'iters', outm.iters, ...
            'errQLen', errQ, 'speedup', t_exact / t_mcd, ...
            'speedup_e2e', (t_reach + t_exact) / (t_reach + t_mcd));
        rows(end + 1) = r; %#ok<AGROW>

        if options.verbose
            line_printf('%3d %4d %9d %9d %8d %9d %9.3f %9.3f %9.3f %6d %10.2e %9.2f %9.2f\n', ...
                K, N, nS, r.sumMk, r.maxMk, nnzQ, t_reach, t_exact, t_mcd, ...
                outm.iters, errQ, r.speedup, r.speedup_e2e);
        end
    end
end

T = struct2table(rows);
end

% ------------------------------------------------------------------------
function n = i_numstates(N, K)
% compositions of N over K stations
n = round(exp(gammaln(N + K) - gammaln(N + 1) - gammaln(K)));
end

% ------------------------------------------------------------------------
function [mu, P, servers] = i_case(K)
% single-class closed cyclic network: one delay (think) station, K-1 queues
mu = 2 + mod(0:(K - 1), 3);          % deterministic, well-conditioned rates
mu(1) = 1;
servers = ones(1, K); servers(1) = Inf;
P = zeros(K);
for i = 1:K, P(i, mod(i, K) + 1) = 1; end
end

% ------------------------------------------------------------------------
function T = i_next(s, ii, jj)
E = numel(ii);
T = zeros(E, numel(s));
m = 0;
for a = 1:E
    i = ii(a);
    if s(i) > 0
        t = s; t(i) = t(i) - 1; t(jj(a)) = t(jj(a)) + 1;
        m = m + 1;
        T(m, :) = t;
    end
end
T = T(1:m, :);
end

% ------------------------------------------------------------------------
function Q = i_solverctmc(mu, servers, N)
K = numel(mu);
model = Network('bench');
st = cell(1, K);
st{1} = Delay(model, 'Think');
for i = 2:K
    st{i} = Queue(model, sprintf('Q%d', i), SchedStrategy.PS);
    st{i}.setNumberOfServers(servers(i));
end
job = ClosedClass(model, 'Jobs', N, st{1});
for i = 1:K, st{i}.setService(job, Exp(mu(i))); end
model.link(Network.serialRouting(st{:}));
Q = SolverCTMC(model).avgTable();
end
