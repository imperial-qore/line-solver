function [QN, UN, RN, TN, CN, XN, totiter] = solver_mam_ldqbd(sn, options)
% SOLVER_MAM_LDQBD Solve single-class closed queueing networks using LD-QBD
%
% Uses Level-Dependent Quasi-Birth-Death (LD-QBD) process to compute
% exact performance metrics for single-class closed queueing networks
% consisting of a Delay (infinite server) and a Queue (FCFS).
%
% The LD-QBD approach models the system where:
%   - Level n = number of jobs at the queue (0 <= n <= N)
%   - Jobs at delay = N - n
%   - Transition rates depend on the current level
%
% Supports PH-type service distributions (Exp, Erlang, HyperExp, etc.)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Validate model structure
M = sn.nstations;
K = sn.nclasses;
N = sn.njobs';

% Check: single-class closed network
if K ~= 1
    line_error(mfilename, 'LDQBD method requires a single-class model.');
end

if ~isfinite(N)
    line_error(mfilename, 'LDQBD method requires a closed model.');
end

% Check: must have exactly one delay and one queue
nDelay = sum(sn.sched == SchedStrategy.INF);
nQueue = sum(sn.sched == SchedStrategy.FCFS);

if nDelay ~= 1 || nQueue ~= 1 || M ~= 2
    line_error(mfilename, 'LDQBD method requires exactly one Delay and one Queue station.');
end

%% Identify stations
delayIdx = find(sn.sched == SchedStrategy.INF);
queueIdx = find(sn.sched == SchedStrategy.FCFS);

%% Get service parameters
PH = sn.proc;
rates = sn.rates;
nservers = sn.nservers;

% Delay station rate (infinite server)
lambda_d = rates(delayIdx, 1);

% Get routing probability from Delay to Queue
% sn.rt is (M*K) x (M*K) routing table
rt = sn.rt;
p_delay_to_queue = rt(delayIdx, queueIdx);  % Since K=1, simple indexing works

% Effective arrival rate to queue = delay rate * routing probability
lambda_eff = lambda_d * p_delay_to_queue;

% Queue service process
PH_queue = PH{queueIdx}{1};
nServers = nservers(queueIdx);

% Check if queue service is exponential (1x1 matrix) or PH
if numel(PH_queue{1}) == 1
    % Exponential service
    mu = -PH_queue{1};  % Rate from D0 matrix
    nPhases = 1;
    isPH = false;
else
    % PH-type service (Erlang, HyperExp, etc.)
    nPhases = size(PH_queue{1}, 1);
    isPH = true;
    D0 = PH_queue{1};  % Sub-generator (negative diag = rates out)
    D1 = PH_queue{2};  % Absorption/completion transitions
    alpha = map_pie(PH_queue);  % Initial phase distribution
end

%% Construct LD-QBD matrices
% For a closed network with N jobs:
%   Level n = number of jobs at queue (0 <= n <= N)
%   Q0^(n): upward transitions (arrival from delay)
%   Q1^(n): local transitions (within level)
%   Q2^(n): downward transitions (departure from queue)

Q0 = cell(N, 1);   % {Q0^(0), Q0^(1), ..., Q0^(N-1)}
Q1 = cell(N+1, 1); % {Q1^(0), Q1^(1), ..., Q1^(N)}
Q2 = cell(N, 1);   % {Q2^(1), Q2^(2), ..., Q2^(N)}

if ~isPH
    % Exponential service case (scalar matrices)
    for n = 0:N-1
        % Arrival rate: (N-n) * lambda_eff (effective rate to queue)
        Q0{n+1} = (N - n) * lambda_eff;
    end

    for n = 0:N
        arrival_rate = (N - n) * lambda_eff;
        if n > 0
            % Service rate: min(n, nServers) * mu
            departure_rate = min(n, nServers) * mu;
        else
            departure_rate = 0;
        end
        Q1{n+1} = -(arrival_rate + departure_rate);
    end

    for n = 1:N
        % Service rate: min(n, nServers) * mu
        Q2{n} = min(n, nServers) * mu;
    end
else
    % PH-type service case (matrix-valued)
    % State space at level n (n>0): nPhases phases
    % State space at level 0: 1 state (empty queue)

    % Level 0 -> 1: arrival starts service in some phase
    % 1 x nPhases matrix: arrival rate * initial phase distribution
    Q0{1} = N * lambda_eff * alpha;

    % Level n -> n+1 (n >= 1): arrival, preserve phase
    for n = 1:N-1
        % nPhases x nPhases matrix: diagonal (preserve phase)
        Q0{n+1} = (N - n) * lambda_eff * eye(nPhases);
    end

    % Level 0: 1x1 (no service, only arrivals)
    Q1{1} = -(N * lambda_eff);

    % Level n >= 1: phase transitions within level
    for n = 1:N
        arrival_rate = (N - n) * lambda_eff;
        % Number of active servers at level n
        c_n = min(n, nServers);

        % Local transition matrix: D0 (phase transitions) - arrivals
        % Scale D0 by number of active servers for multi-server
        Q1{n+1} = c_n * D0 - arrival_rate * eye(nPhases);
    end

    % Level 1 -> 0: service completion, go to empty state
    % nPhases x 1 column vector (sum over destination phases)
    c_1 = min(1, nServers);
    Q2{1} = c_1 * D1 * ones(nPhases, 1);

    % Level n -> n-1 (n >= 2): service completion, next job starts
    % D1(i,j) is rate of completing from phase i and next starting in phase j
    for n = 2:N
        c_n = min(n, nServers);
        % nPhases x nPhases: D1 encodes completion and restart phase
        Q2{n} = c_n * D1;
    end
end

%% Solve LD-QBD using the ldqbd function
ldqbd_options = struct('epsilon', options.tol, 'maxIter', options.iter_max, 'verbose', false);

[R, pi_ldqbd] = ldqbd(Q0, Q1, Q2, ldqbd_options);

%% Compute performance metrics from steady-state distribution
% Mean queue length at the queue station
mean_queue = (0:N) * pi_ldqbd';

% Mean queue length at delay station (complementary)
mean_delay = N - mean_queue;

%% Compute throughput and utilization
% For delay (infinite server): each job at delay generates arrivals to queue
% at the effective rate lambda_eff = lambda_d * P(Delay->Queue)
% System throughput X = mean_delay * lambda_eff (flow balance)
X = mean_delay * lambda_eff;

% Mean service time at queue
if ~isPH
    mean_service = 1/mu;
else
    mean_service = map_mean(PH_queue);
end

% Queue utilization: probability queue is busy (at least 1 job)
% For multi-server: average fraction of server capacity in use
if nServers == 1
    util_queue = 1 - pi_ldqbd(1);
else
    % For c-server queue: U = sum_{n=1}^{N} min(n,c)/c * pi(n)
    util_queue = 0;
    for n = 1:N
        util_queue = util_queue + (min(n, nServers) / nServers) * pi_ldqbd(n+1);
    end
end

% Response time at queue (using Little's law: R = Q/X)
if X > 0
    R_queue = mean_queue / X;
else
    R_queue = 0;
end

% Response time at delay (constant for infinite server)
R_delay = 1 / lambda_d;

%% Populate output matrices
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

% Delay station metrics
QN(delayIdx, 1) = mean_delay;
UN(delayIdx, 1) = mean_delay;  % For infinite server, U = Q
RN(delayIdx, 1) = R_delay;
TN(delayIdx, 1) = X;

% Queue station metrics
QN(queueIdx, 1) = mean_queue;
UN(queueIdx, 1) = util_queue / nServers;  % Per-server utilization
RN(queueIdx, 1) = R_queue;
TN(queueIdx, 1) = X;

% System-level metrics
XN(1) = X;
CN(1) = R_delay + R_queue;  % Cycle time

totiter = 1;  % LDQBD is a direct method

end
