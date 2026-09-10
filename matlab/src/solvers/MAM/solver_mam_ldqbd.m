function [QN, UN, RN, TN, CN, XN, totiter, ld] = solver_mam_ldqbd(sn, options)
% SOLVER_MAM_LDQBD Solve single-class Delay/Queue networks using LD-QBD
%
% Uses a Level-Dependent Quasi-Birth-Death (LD-QBD) process to compute
% performance metrics for single-class networks with one infinite-server station
% and one FCFS Queue (multi-server, PH service supported).
%
% Exactness: exact for exponential service at any number of servers, and for PH
% service at any number of servers. The multiserver PH chain is built by
% LDQBD_MPHC, whose level coordinate is the MULTISET of the phases the min(n,c)
% busy servers sit in; the collapsed single-phase approximation this solver used
% until 2026-08-18 (one PH process run at min(n,c) times its speed, ~1e-2
% relative against SolverCTMC) is gone.
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
% The shape rules live in MAM_LDQBD_APPLICABLE, which SolverMAM.supportsModelMethod
% and solver_mam_analyzer ask as well: one body, so the report and this run agree.
[shapeOk, shapeWhy, regime] = mam_ldqbd_applicable(sn);
if ~shapeOk
    line_error(mfilename, shapeWhy);
end
M = sn.nstations;
K = sn.nclasses;
N = sn.njobs';

isOpen = strcmp(regime, 'open');

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

%% Setup and delay-off at the queue
% The station alternates OFF -> setup -> busy -> delay-off around the service,
% so the chain carries phases the block builder below has no place for. The
% closed regime hands the whole chain to qbd_setupdelayoff_closed; every other
% combination is refused BY NAME rather than answered as if the server were
% always warm, which is what this solver did until 2026-09 and is BUG-78.
hasSetup = isfield(sn,'hassetup') && numel(sn.hassetup) >= queueIdx && sn.hassetup(queueIdx);
alpharate = NaN; alphascv = NaN; betarate = NaN; betascv = NaN;
if hasSetup
    if isOpen
        line_error(mfilename, ['Open LDQBD does not model a setup/delay-off server; ' ...
            'use options.method=''dec.source'', whose qbd_setupdelayoff covers the open case.']);
    end
    if isPH || nservers(queueIdx) > 1
        line_error(mfilename, ['Closed LDQBD models a setup/delay-off server with ' ...
            'exponential service at a single server only; this station has ' ...
            'phase-type service or several servers.']);
    end
    nodeIdx = sn.stationToNode(queueIdx);
    np = sn.nodeparam{nodeIdx};
    fparam = np{end};
    alpharate = map_lambda(fparam.setupTime);
    alphascv = map_scv(fparam.setupTime);
    betarate = map_lambda(fparam.delayoffTime);
    betascv = map_scv(fparam.delayoffTime);
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
    % Peak capacity normalizes the utilization, and it is the LARGEST factor the
    % table declares, not the saturated one: a non-monotone alpha peaks in the
    % middle. Same rule as CTMC's ceff = max(nservers, max(lldscaling(ist,:))),
    % which is what makes the two report the same number.
    utilPeak = max(nServers, max(lld));
else
    lld = [];
    lldlimit = 0;
    sfMax = nServers;
    utilPeak = nServers;
end
if hasSetup && hasLLD
    line_error(mfilename, ['Closed LDQBD models a setup/delay-off server at its nominal ' ...
        'rate only; this station also declares a load-dependent scaling.']);
end

%% Arrival rate per level and number of levels
rt = sn.rt;
if isOpen
    % External Poisson arrivals only; MAM_LDQBD_APPLICABLE refused a MAP above.
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
    % PH service: the level carries the MULTISET of the phases the min(n,c)
    % busy servers sit in, which is exact at any number of servers. At c = 1
    % the multiset is just the phase, so this reproduces the single-server
    % blocks (sf(n)*D0 - arr*I, sf(n)*D1) entry for entry.
    [Q0, Q1, Q2] = ldqbd_mphc(D0, D1, alpha, nServers, arrRate, sf);
end

%% Solve LD-QBD
if hasSetup
    % SETUP AND DELAY-OFF, the closed vacation queue. The level-dependent chain
    % this needs is the one above with two extra phase families -- the setup
    % above level 0 and the delay-off at level 0 -- and qbd_setupdelayoff_closed
    % builds and solves exactly that, so it is called rather than duplicated.
    % Without it the blocks above describe a server that is ALWAYS warm and the
    % answer is byte-identical across any setup mean (BUG-78).
    [mean_queue, X_setup] = qbd_setupdelayoff_closed(N, 1/lambda_eff, mu, ...
        alpharate, alphascv, betarate, betascv);
    pi_ldqbd = [];
else
    ldqbd_options = struct('epsilon', options.tol, 'maxIter', options.iter_max, 'verbose', false);
    [R, pi_ldqbd] = ldqbd(Q0, Q1, Q2, ldqbd_options); %#ok<ASGLU>
    %% Performance metrics (per-level stationary distribution pi_ldqbd)
    mean_queue = (0:Nlev) * pi_ldqbd';
end

% Utilization is the fraction of the station's PEAK capacity in use,
% sum_n pi(n) * sf(n) / utilPeak, which is the work-based convention CTMC, MVA,
% NC and serial SSA all report. utilPeak = max(c, max(alpha)) is CTMC's own
% normalizer, so the two agree state for state.
%
% This subsumes the two special cases it replaces rather than approximating
% them: without load dependence sf(n) = min(n,c) and utilPeak = c, giving the
% average fraction of c servers in use; at c = 1 that is sf(n) = 1 for every
% n >= 1, so the sum collapses to 1 - pi(0).
%
% It used to report 1 - pi(0) under load dependence, i.e. P(busy). That reads a
% station running alpha(n) times faster as no busier than one running at its
% nominal rate, and put MAM 0.9587 against CTMC's 0.6612 on a 4-job closed model
% with alpha = [1 1.5 2 2.5].
util_ps = 0;
if hasSetup
    % With a setup the server is DELIVERING work only in the busy phase, so the
    % level occupancy over-counts it: a level is occupied during the setup too.
    % The utilization law gives the same work-based number without needing the
    % per-phase vector, X*E[S]/peak, which is what the level sum reduces to
    % without a setup.
    util_ps = X_setup * mean_service / utilPeak;
else
    for n = 1:Nlev
        util_ps = util_ps + (sf(n) / utilPeak) * pi_ldqbd(n+1);
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
        'sf', {sf}, 'utilPeak', utilPeak, ...
        'lambda_eff', lambda_eff, 'delayRate', delayRate, 'N', Npop);
end

end
