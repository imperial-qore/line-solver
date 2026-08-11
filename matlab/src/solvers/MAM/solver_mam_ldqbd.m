function [QN, UN, RN, TN, CN, XN, totiter, ld] = solver_mam_ldqbd(sn, options)
% SOLVER_MAM_LDQBD Solve single-class Delay/Queue networks using LD-QBD
%
% Uses a Level-Dependent Quasi-Birth-Death (LD-QBD) process to compute
% performance metrics for single-class networks with one infinite-server station
% and one FCFS Queue (multi-server, PH service supported).
%
% Exactness: exact for exponential service at any number of servers, and for PH
% service at a single server. For PH service with c > 1 servers it is an
% approximation: the c parallel PH servers are collapsed into one PH process
% scaled by min(n,c), which ignores the phase of each individual busy server
% (the exact chain tracks the multiset of the min(n,c) in-service phases).
% Measured against SolverCTMC the residual error is ~1e-2 relative on Erlang and
% HyperExp, still ~16x below the dec.source error on the same models.
%
% Two regimes are handled:
%   CLOSED: one Delay (INF) + one Queue, finite population N. Level n = jobs at
%           the queue (0 <= n <= N); arrival rate from the delay is (N-n)*lambda.
%   OPEN:   one Source (EXT) + one Queue, open class (Poisson arrivals). Level
%           n = jobs at the queue, truncated at M_trunc; arrival rate is the
%           constant external rate lambda. M_trunc is taken from options.cutoff
%           or chosen so the truncated tail probability is negligible.
%
% Both regimes share the same block-tridiagonal generator, differing only in the
% per-level arrival rate and the top level. Service is min(n,c)*mu (exact M/M/c
% boundary) or its PH generalisation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Validate model structure
M = sn.nstations;
K = sn.nclasses;
N = sn.njobs';

if K ~= 1
    line_error(mfilename, 'LDQBD method requires a single-class model.');
end

nDelay  = sum(sn.sched == SchedStrategy.INF);
nQueue  = sum(sn.sched == SchedStrategy.FCFS);
nSource = sum(sn.sched == SchedStrategy.EXT);

isOpen = ~isfinite(N);
if isOpen
    if nSource ~= 1 || nQueue ~= 1 || M ~= 2
        line_error(mfilename, 'Open LDQBD method requires exactly one Source and one Queue station.');
    end
else
    if nDelay ~= 1 || nQueue ~= 1 || M ~= 2
        line_error(mfilename, 'Closed LDQBD method requires exactly one Delay and one Queue station.');
    end
end

%% Identify stations
queueIdx = find(sn.sched == SchedStrategy.FCFS);
if isOpen
    srcIdx = find(sn.sched == SchedStrategy.EXT);
else
    delayIdx = find(sn.sched == SchedStrategy.INF);
end

%% Service process at the queue
PH = sn.proc;
rates = sn.rates;
nservers = sn.nservers;
PH_queue = PH{queueIdx}{1};
nServers = nservers(queueIdx);

if numel(PH_queue{1}) == 1
    mu = -PH_queue{1};   % exponential service rate
    nPhases = 1;
    isPH = false;
    mean_service = 1/mu;
else
    nPhases = size(PH_queue{1}, 1);
    isPH = true;
    D0 = PH_queue{1};
    D1 = PH_queue{2};
    alpha = map_pie(PH_queue);
    mean_service = map_mean(PH_queue);
end

%% Per-level service factor: load-dependent scaling if set, else min(n,c)
% sn.lldscaling(queueIdx, n) is the load-dependent multiplier on the base rate;
% when absent, a c-server queue scales as min(n,c). This is the level-dependent
% factor that the LD-QBD applies in each downward (departure) block.
hasLLD = sn_has_load_dependence(sn) && ~isempty(sn.lldscaling) ...
    && size(sn.lldscaling, 1) >= queueIdx && any(sn.lldscaling(queueIdx, :) ~= 1);
if hasLLD
    lld = sn.lldscaling(queueIdx, :);
    lldlimit = numel(lld);
    sfMax = lld(lldlimit);          % saturated factor (capacity ceiling)
else
    lld = [];
    lldlimit = 0;
    sfMax = nServers;
end

%% Arrival rate per level and number of levels
rt = sn.rt;
if isOpen
    % External Poisson arrivals only (MAP/MMPP arrivals are not yet supported).
    arrProc = PH{srcIdx}{1};
    if numel(arrProc{1}) > 1
        line_error(mfilename, ['Open LDQBD method currently supports Poisson (exponential) ' ...
            'arrivals only; the Source uses a MAP/MMPP process.']);
    end
    lambda = rates(srcIdx, 1);
    lambda_eff = lambda * rt(srcIdx, queueIdx);
    rho = lambda_eff * mean_service / sfMax;   % sfMax = saturated capacity factor
    if rho >= 1
        line_error(mfilename, sprintf(['Open LDQBD method requires a stable queue ' ...
            '(rho = %.4f >= 1). Increase service capacity or reduce the arrival rate.'], rho));
    end
    % Truncation level: explicit cutoff, else enough levels for a negligible tail.
    if isfield(options, 'cutoff') && isscalar(options.cutoff) && isfinite(options.cutoff)
        Nlev = max(nServers + 1, round(options.cutoff));
    else
        tailTol = 1e-10;
        Nlev = nServers + ceil(log(tailTol) / log(rho));
        Nlev = min(max(Nlev, nServers + 10), 100000);
    end
    arrRate = lambda_eff * ones(1, Nlev + 1); % arrRate(n+1) is the rate out of level n
    arrRate(Nlev + 1) = 0;                    % truncation: no arrivals above the top level
else
    lambda_d = rates(delayIdx, 1);
    lambda_eff = lambda_d * rt(delayIdx, queueIdx);
    Nlev = N;
    arrRate = (N - (0:N)) * lambda_eff;       % finite-source rate (N-n)*lambda_eff, 0 at n=N
end

% Per-level service factor sf(n), n = 1..Nlev
sf = zeros(1, Nlev);
for n = 1:Nlev
    if hasLLD
        sf(n) = lld(min(n, lldlimit));
    else
        sf(n) = min(n, nServers);
    end
end

%% Construct LD-QBD block-tridiagonal generator
%   Q0^(n): upward (arrival)    Q1^(n): local    Q2^(n): downward (departure)
Q0 = cell(Nlev, 1);
Q1 = cell(Nlev + 1, 1);
Q2 = cell(Nlev, 1);

if ~isPH
    for n = 0:Nlev-1
        Q0{n+1} = arrRate(n+1);
    end
    for n = 0:Nlev
        departure_rate = 0;
        if n > 0
            departure_rate = sf(n) * mu;
        end
        Q1{n+1} = -(arrRate(n+1) + departure_rate);
    end
    for n = 1:Nlev
        Q2{n} = sf(n) * mu;
    end
else
    Q0{1} = arrRate(1) * alpha;                 % level 0 -> 1: start service in a phase
    for n = 1:Nlev-1
        Q0{n+1} = arrRate(n+1) * eye(nPhases);  % level n -> n+1: preserve phase
    end
    Q1{1} = -arrRate(1);                        % level 0: only arrivals
    for n = 1:Nlev
        Q1{n+1} = sf(n) * D0 - arrRate(n+1) * eye(nPhases);
    end
    Q2{1} = sf(1) * D1 * ones(nPhases, 1);      % level 1 -> 0: empty the queue
    for n = 2:Nlev
        Q2{n} = sf(n) * D1;                     % level n -> n-1: complete and restart
    end
end

%% Solve LD-QBD
ldqbd_options = struct('epsilon', options.tol, 'maxIter', options.iter_max, 'verbose', false);
[R, pi_ldqbd] = ldqbd(Q0, Q1, Q2, ldqbd_options); %#ok<ASGLU>

%% Performance metrics (per-level stationary distribution pi_ldqbd)
mean_queue = (0:Nlev) * pi_ldqbd';

% Utilization: probability busy for a (load-dependent) single server, else the
% average fraction of c servers in use.
if hasLLD || nServers == 1
    util_ps = 1 - pi_ldqbd(1);
else
    util_ps = 0;
    for n = 1:Nlev
        util_ps = util_ps + (min(n, nServers) / nServers) * pi_ldqbd(n+1);
    end
end

QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

if isOpen
    % Accepted/served throughput = arrival rate minus truncation blocking
    % (= lambda_eff when the truncation tail is negligible).
    X = lambda_eff * (1 - pi_ldqbd(Nlev + 1));
    if X > 0
        R_queue = mean_queue / X;
    else
        R_queue = 0;
    end
    % Source station: pass-through, no queueing.
    QN(srcIdx, 1) = 0;
    UN(srcIdx, 1) = 0;
    RN(srcIdx, 1) = 0;
    TN(srcIdx, 1) = X;
    % Queue station.
    QN(queueIdx, 1) = mean_queue;
    UN(queueIdx, 1) = util_ps;
    RN(queueIdx, 1) = R_queue;
    TN(queueIdx, 1) = X;
    XN(1) = X;
    CN(1) = R_queue;
else
    mean_delay = N - mean_queue;
    X = mean_delay * lambda_eff;
    if X > 0
        R_queue = mean_queue / X;
    else
        R_queue = 0;
    end
    R_delay = 1 / rates(delayIdx, 1);
    % Delay throughput is not the queue flow when rt<1; see _kb/06-solver-catalog.md for rationale
    QN(delayIdx, 1) = mean_delay;
    UN(delayIdx, 1) = mean_delay;       % infinite server: U = Q
    RN(delayIdx, 1) = R_delay;
    TN(delayIdx, 1) = mean_delay * lambda_d;
    % Queue station metrics
    QN(queueIdx, 1) = mean_queue;
    UN(queueIdx, 1) = util_ps;
    RN(queueIdx, 1) = R_queue;
    TN(queueIdx, 1) = X;
    XN(1) = X;
    CN(1) = R_delay + R_queue;
end

totiter = 1;  % LDQBD is a direct method

%% Optional: expose the LD-QBD blocks and parameters (for the SolverENV
% state-vector analyzer's MAM backend). Built only when requested.
if nargout >= 8
    if isOpen
        refIdx = srcIdx;
        delayRate = NaN;
        Npop = Inf;
    else
        refIdx = delayIdx;
        delayRate = rates(delayIdx, 1);
        Npop = N;
    end
    ld = struct('Q0', {Q0}, 'Q1', {Q1}, 'Q2', {Q2}, ...
        'Nlev', Nlev, 'nPhases', nPhases, 'isPH', isPH, 'isOpen', isOpen, ...
        'queueIdx', queueIdx, 'refIdx', refIdx, 'M', M, ...
        'nServers', nServers, 'mean_service', mean_service, 'hasLLD', hasLLD, ...
        'lambda_eff', lambda_eff, 'delayRate', delayRate, 'N', Npop);
end

end
