function [QN, UN, RN, TN, CN, XN, iter] = solver_mam_ag(sn, options)
% SOLVER_MAM_AG AG methods for SolverMAM
%
% [QN, UN, RN, TN, CN, XN, ITER] = SOLVER_MAM_AG(SN, OPTIONS)
%
% Uses RCAT (Reversed Compound Agent Theorem) to find product-form
% solutions for queueing networks.
%
% Methods:
%   'inap'     - Iterative Numerical Approximation Procedure (default, fast)
%   'inapplus' - Improved INAP with weighted rates (no normalization)
%   'exact'    - Not available (autocat moved to line-legacy.git)
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;

% Set default max states for truncation
if isfield(options, 'config') && isfield(options.config, 'maxStates')
    maxStates = options.config.maxStates;
else
    maxStates = 100;
end

% Set default tolerances
if isfield(options, 'iter_tol') && ~isempty(options.iter_tol)
    tol = options.iter_tol;
else
    tol = 1e-6;
end

if isfield(options, 'iter_max') && ~isempty(options.iter_max)
    maxiter = options.iter_max;
else
    maxiter = 1000;
end

% Build RCAT model from network structure
[R, AP, processMap, actionMap, N] = build_rcat(sn, maxStates);

% Check if we have a valid model
numProcesses = max(processMap(:));
numActions = size(AP, 1);

% Return early only if no processes found
if numProcesses == 0
    line_warning(mfilename, 'Network could not be mapped to RCAT format (no processes found).\n');
    QN = zeros(M, K);
    UN = zeros(M, K);
    RN = zeros(M, K);
    TN = zeros(M, K);
    CN = zeros(1, K);
    XN = zeros(1, K);
    iter = 0;
    return;
end

% If no actions but we have processes, solve using local rates only
% This handles single-queue G-networks (Source -> Queue -> Sink)
if numActions == 0
    % No inter-station actions: solve equilibrium using only L matrices
    x = [];
    pi = cell(1, numProcesses);
    Q = cell(1, numProcesses);
    for p = 1:numProcesses
        L = R{1, p};  % Local rate matrix is in R{numActions+1, p} = R{1, p} when numActions=0
        % Convert to valid generator matrix
        Qp = L - diag(L * ones(size(L, 1), 1));
        Q{p} = ctmc_makeinfgen(Qp);
        % Solve for equilibrium
        pi{p} = ctmc_solve(Q{p});
    end
    iter = 0;
    [QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N);
    return;
end

% Choose solver method
method = options.method;
if strcmp(method, 'default')
    method = 'inap';
end

switch method
    case 'inap'
        % Fast iterative heuristic
        [x, pi, Q, iter] = inap(R, AP, tol, maxiter, 'inap');

    case 'inapplus'
        % Improved INAP with weighted rates (no normalization)
        [x, pi, Q, iter] = inap(R, AP, tol, maxiter, 'inapplus');

    case 'exact'
        % Optimization-based solver using autocat (not available in this version)
        line_warning(mfilename, '''exact'' method not available. Falling back to inap.\n');
        [x, pi, Q, iter] = inap(R, AP, tol, maxiter, 'inap');

    otherwise
        line_error(mfilename, 'Unknown method: %s\n', method);
end

% Convert RCAT solution to LINE metrics
[QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N);

end

%% Local Functions

function [x, pi, Q, iter] = inap(R, AP, tol, maxiter, method)
% INAP Iterative Numerical Approximation Procedure for RCAT
%
% Methods:
%   'inap':     x(a) = mean(Aa(i,j) * pi(i) / pi(j))
%   'inapplus': x(a) = sum(Aa(i,j) * pi(i))

if nargin < 5 || isempty(method)
    method = 'inap';
end

% Parse R and AP
A = size(AP, 1);
ACT = AP(:, 1);
PSV = AP(:, 2);
numProcesses = max(AP(:));

% Extract rate matrices
Aa = cell(1, A);
Pb = cell(1, A);
for a = 1:A
    Aa{a} = R{a, 1};
    Pb{a} = R{a, 2};
end

% Extract local rates
L = cell(1, numProcesses);
for k = 1:numProcesses
    L{k} = R{A+1, k};
end

% Get state space sizes
N = zeros(1, numProcesses);
for k = 1:numProcesses
    N(k) = size(L{k}, 1);
end

% Initialize action rates randomly
x = rand(A, 1);

% Compute initial equilibrium
[pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N);

iter = 1;
while iter <= maxiter
    % Save old iteration result
    piprev = pi;

    % Update each action rate
    for a = 1:A
        k = ACT(a);

        if strcmp(method, 'inapplus')
            % inapplus: LAMBDA(i,j) = Aa{a}(i,j) * pi{k}(i)
            % x(a) = sum(LAMBDA) for non-zero entries
            LAMBDA_sum = 0;
            for i = 1:N(k)
                for j = 1:N(k)
                    if Aa{a}(i,j) > 0
                        LAMBDA_sum = LAMBDA_sum + Aa{a}(i,j) * pi{k}(i);
                    end
                end
            end
            if LAMBDA_sum > 0
                x(a) = LAMBDA_sum;
            end
        else
            % inap: LAMBDA(i,j) = Aa{a}(i,j) * pi{k}(i) / pi{k}(j)
            % x(a) = mean(LAMBDA) for non-zero entries
            LAMBDA_vec = [];
            for i = 1:N(k)
                for j = 1:N(k)
                    if Aa{a}(i,j) > 0 && pi{k}(j) > 0
                        LAMBDA_vec(end+1) = Aa{a}(i,j) * pi{k}(i) / pi{k}(j);
                    end
                end
            end
            if ~isempty(LAMBDA_vec)
                x(a) = mean(LAMBDA_vec);
            end
        end
    end

    % Recompute equilibrium with new x
    [pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N);

    % Convergence check
    maxErr = 0;
    for k = 1:numProcesses
        maxErr = max([maxErr, norm(pi{k} - piprev{k}, 1)]);
    end

    if maxErr < tol
        return
    end

    iter = iter + 1;
end

end

function [pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N)
% Compute equilibrium distribution for each process given action rates x

Q = cell(1, numProcesses);
pi = cell(1, numProcesses);

for k = 1:numProcesses
    % Start with local/hidden rates
    Qk = L{k} - diag(L{k} * ones(N(k), 1));

    % Add contributions from each action
    for c = 1:A
        if PSV(c) == k
            % Process k is passive for action c: add x(c) * Pb{c}
            Qk = Qk + x(c) * Pb{c} - diag(Pb{c} * ones(N(k), 1));
        elseif ACT(c) == k
            % Process k is active for action c: add Aa{c}
            Qk = Qk + Aa{c} - diag(Aa{c} * ones(N(k), 1));
        end
    end

    % Convert to valid infinitesimal generator
    Q{k} = ctmc_makeinfgen(Qk);

    % Solve for equilibrium distribution using birth-death recursion for
    % tridiagonal generators (numerically stable for large state spaces),
    % falling back to ctmc_solve otherwise
    if is_tridiagonal(Q{k})
        pi{k} = birth_death_solve(Q{k});
    else
        pi{k} = ctmc_solve(Q{k});
    end
end

end

function [R, AP, processMap, actionMap, N] = build_rcat(sn, maxStates)
% BUILD_RCAT Convert LINE network structure to RCAT format

if nargin < 2
    maxStates = 100;
end

M = sn.nstations;
K = sn.nclasses;
rt = sn.rt;  % (M*K) x (M*K) routing table

% Identify station types
sourceStations = [];
sinkStations = [];
queueStations = [];

for ist = 1:M
    nodeIdx = sn.stationToNode(ist);
    if sn.nodetype(nodeIdx) == NodeType.Source
        sourceStations(end+1) = ist;
    elseif sn.nodetype(nodeIdx) == NodeType.Sink
        sinkStations(end+1) = ist;
    else
        % Queue, Delay, or other service stations
        queueStations(end+1) = ist;
    end
end

% Create process mapping: each (station, class) pair at queue stations
% Note: Signal classes (negative customers) don't create separate processes
% as they only modify the state of positive customer processes
processIdx = 0;
processMap = zeros(M, K);
for ist = queueStations
    for r = 1:K
        % Skip Signal classes - they don't have their own queue state
        if sn.issignal(r)
            continue;
        end
        % Check if this station serves this class
        if ~isnan(sn.rates(ist, r)) && sn.rates(ist, r) > 0
            processIdx = processIdx + 1;
            processMap(ist, r) = processIdx;
        end
    end
end
numProcesses = processIdx;

if numProcesses == 0
    R = {};
    AP = [];
    actionMap = [];
    N = [];
    return;
end

% Check for signal classes - G-networks with signals
% G-network signals (negative customers, catastrophes, batch removal) are handled by:
% 1. Adding signal arrival effects to the L matrix for positive customer processes
% 2. Creating actions for signal routing between stations (job removal at destination)

% Determine number of states for each process
N = zeros(1, numProcesses);
for p = 1:numProcesses
    [ist, r] = find(processMap == p);
    if ~isempty(ist)
        ist = ist(1); r = r(1);
        if sn.njobs(r) < Inf  % Closed class
            N(p) = sn.njobs(r) + 1;  % States 0, 1, ..., njobs
        else  % Open class
            N(p) = maxStates;  % Truncate at maxStates
        end
    end
end

% Count actions: each routing transition (i,r) -> (j,s) where P > 0
actionIdx = 0;
actionMap = struct('from_station', {}, 'from_class', {}, ...
                   'to_station', {}, 'to_class', {}, 'prob', {}, ...
                   'isNegative', {}, 'isCatastrophe', {}, 'removalDistribution', {});

for ist = queueStations
    for r = 1:K
        if processMap(ist, r) > 0
            % Check if class r is a negative signal class
            isNegativeClass = false;
            isCatastropheClass = false;
            removalDist = [];
            if sn.issignal(r) && ~isnan(sn.signaltype{r}) && sn.signaltype{r} == SignalType.NEGATIVE
                isNegativeClass = true;
                % Check if class r is a catastrophe signal
                if isfield(sn, 'isCatastrophe') && ~isempty(sn.isCatastrophe) && sn.isCatastrophe(r) > 0
                    isCatastropheClass = true;
                end
                % Get removal distribution for this class
                if isfield(sn, 'signalRemovalDist') && ~isempty(sn.signalRemovalDist) && r <= length(sn.signalRemovalDist)
                    removalDist = sn.signalRemovalDist{r};
                end
            end

            for jst = queueStations
                for s = 1:K
                    if processMap(jst, s) > 0
                        % Get routing probability
                        prob_ij_rs = rt((ist-1)*K + r, (jst-1)*K + s);
                        if prob_ij_rs > 0 && (ist ~= jst || r ~= s)
                            % This is an action (departure from i,r triggers arrival at j,s)
                            actionIdx = actionIdx + 1;
                            actionMap(actionIdx).from_station = ist;
                            actionMap(actionIdx).from_class = r;
                            actionMap(actionIdx).to_station = jst;
                            actionMap(actionIdx).to_class = s;
                            actionMap(actionIdx).prob = prob_ij_rs;
                            actionMap(actionIdx).isNegative = isNegativeClass;
                            actionMap(actionIdx).isCatastrophe = isCatastropheClass;
                            actionMap(actionIdx).removalDistribution = removalDist;
                        end
                    end
                end
            end
        end
    end
end
numActions = actionIdx;

% Initialize R and AP
R = cell(numActions + 1, max(numProcesses, 2));
if numActions > 0
    AP = zeros(numActions, 2);
else
    AP = zeros(0, 2);  % Empty matrix when no actions
end

% Identify sink nodes (nodetype = -1 = NodeType.Sink)
% Use row vector to ensure for-loop doesn't execute when empty
sinkNodes = find(sn.nodetype == NodeType.Sink)';
if isempty(sinkNodes)
    sinkNodes = [];  % Ensure empty row vector, not column
end

% Build local/hidden rate matrices L for each process (R{end,k})
for p = 1:numProcesses
    [ist, r] = find(processMap == p);
    ist = ist(1); r = r(1);
    R{numActions + 1, p} = build_local_rates(sn, ist, r, N(p), rt, sourceStations, sinkNodes, K);
end

% Build active and passive matrices for each action
for a = 1:numActions
    am = actionMap(a);

    % Active process (departure)
    ist = am.from_station;
    r = am.from_class;
    p_active = processMap(ist, r);
    AP(a, 1) = p_active;

    % Get service rate
    mu_ir = sn.rates(ist, r);
    prob = am.prob;

    % Active matrix: transition n -> n-1 with rate mu*prob (service completion)
    Aa = zeros(N(p_active));
    for n = 2:N(p_active)
        Aa(n, n-1) = mu_ir * prob;
    end
    % Boundary: at max capacity
    Aa(N(p_active), N(p_active)) = mu_ir * prob;
    R{a, 1} = Aa;

    % Passive process (arrival or signal effect)
    jst = am.to_station;
    s = am.to_class;
    p_passive = processMap(jst, s);
    AP(a, 2) = p_passive;

    Pb = zeros(N(p_passive));
    if am.isNegative
        % NEGATIVE: Job removal at destination (G-network negative customer)
        if am.isCatastrophe
            % CATASTROPHE: All jobs are removed - all states transition to 1 (empty)
            for n = 1:N(p_passive)
                Pb(n, 1) = 1;
            end
        elseif ~isempty(am.removalDistribution)
            % BATCH REMOVAL: Remove a random number of jobs based on distribution
            % P[n, m] = probability of transition from n to m jobs
            dist = am.removalDistribution;
            for n = 1:N(p_passive)
                if n == 1
                    % Empty queue: no effect
                    Pb(1, 1) = 1;
                else
                    % For each possible resulting state m (from 1 to n)
                    for m = 1:n
                        k = n - m;  % Number of jobs to remove to go from n to m
                        % Probability of removing exactly k jobs when queue has n-1 jobs (0-indexed)
                        if m > 1
                            % Remove exactly k jobs: P(removal = k)
                            prob = dist.evalPMF(k);
                            if prob > 0
                                Pb(n, m) = Pb(n, m) + prob;
                            end
                        else
                            % Remove all jobs (m = 1, i.e., state 0): P(removal >= n-1)
                            % = 1 - CDF(n-2) = 1 - sum_{j=0}^{n-2} P(removal = j)
                            cdfNMinus1 = 0;
                            for j = 0:(n-2)
                                cdfNMinus1 = cdfNMinus1 + dist.evalPMF(j);
                            end
                            probAtLeastN = 1 - cdfNMinus1;
                            if probAtLeastN > 0
                                Pb(n, 1) = Pb(n, 1) + probAtLeastN;
                            end
                        end
                    end
                end
            end
        else
            % DEFAULT: Remove exactly 1 job (original behavior)
            % Empty queue: no effect (state 1 stays at state 1)
            Pb(1, 1) = 1;
            % Non-empty queues: decrement (n -> n-1)
            for n = 2:(N(p_passive) - 1)
                Pb(n, n-1) = 1;
            end
            % Boundary at max capacity: decrement
            if N(p_passive) > 1
                Pb(N(p_passive), N(p_passive) - 1) = 1;
            end
        end
    else
        % POSITIVE: Normal job arrival at destination
        for n = 1:(N(p_passive) - 1)
            Pb(n, n+1) = 1;
        end
        % Boundary: at max capacity
        Pb(N(p_passive), N(p_passive)) = 1;
    end
    R{a, 2} = Pb;
end

end

function L = build_local_rates(sn, ist, r, Np, rt, sourceStations, sinkNodes, K)
% Build local/hidden transition matrix for process at station ist, class r
% Note: sinkNodes contains node indices (not station indices) for Sink nodes

L = zeros(Np);

% External arrivals from source - separate positive, negative (single), batch, and catastrophe
lambda_ir_pos = 0;  % Positive arrivals
lambda_ir_neg_single = 0;  % Negative arrivals with single removal (default)
lambda_ir_catastrophe = 0;  % Catastrophe arrivals (remove all)
% Batch removal arrivals: cell array of {rate, distribution} pairs
batchArrivals = {};

for isrc = sourceStations
    for s_src = 1:K
        % Check if source class s_src is a signal
        isSignal = sn.issignal(s_src);

        if isSignal
            % For signals: they route to themselves (Signal -> Signal), but their effect
            % is on positive customers at the destination station. We check if the signal
            % routes to ANY class at this station (not just class r).
            prob_src = 0;
            for s_dst = 1:K
                prob_src = prob_src + rt((isrc-1)*K + s_src, (ist-1)*K + s_dst);
            end
        else
            % For regular classes: direct routing to (ist, r)
            prob_src = rt((isrc-1)*K + s_src, (ist-1)*K + r);
        end

        if prob_src > 0 && ~isnan(sn.rates(isrc, s_src))
            srcRate = sn.rates(isrc, s_src);
            % Check if source class s_src is a negative or catastrophe signal
            if isSignal && ~isnan(sn.signaltype{s_src}) && ...
                    (sn.signaltype{s_src} == SignalType.NEGATIVE || sn.signaltype{s_src} == SignalType.CATASTROPHE)
                % Check if it's a catastrophe (either via isCatastrophe flag or signaltype)
                isCat = (isfield(sn, 'isCatastrophe') && ~isempty(sn.isCatastrophe) && sn.isCatastrophe(s_src)) || ...
                        sn.signaltype{s_src} == SignalType.CATASTROPHE;
                if isCat
                    lambda_ir_catastrophe = lambda_ir_catastrophe + srcRate * prob_src;
                else
                    % Check if it has a removal distribution
                    removalDist = [];
                    if isfield(sn, 'signalRemovalDist') && ~isempty(sn.signalRemovalDist) && s_src <= length(sn.signalRemovalDist)
                        removalDist = sn.signalRemovalDist{s_src};
                    end
                    if ~isempty(removalDist)
                        batchArrivals{end+1} = {srcRate * prob_src, removalDist};
                    else
                        lambda_ir_neg_single = lambda_ir_neg_single + srcRate * prob_src;
                    end
                end
            else
                lambda_ir_pos = lambda_ir_pos + srcRate * prob_src;
            end
        end
    end
end

% Positive arrival transitions: n -> n+1
if lambda_ir_pos > 0
    for n = 1:(Np-1)
        L(n, n+1) = lambda_ir_pos;
    end
end

% Catastrophe arrival transitions: n -> 1 (for all n > 1)
if lambda_ir_catastrophe > 0
    for n = 2:Np
        L(n, 1) = L(n, 1) + lambda_ir_catastrophe;
    end
end

% Batch removal arrival transitions: n -> m at rate λ * P(remove n-m) for m < n
for b = 1:length(batchArrivals)
    rate = batchArrivals{b}{1};
    dist = batchArrivals{b}{2};
    for n = 2:Np
        for m = 1:n
            k = n - m;  % Number of jobs to remove
            if m > 1
                % Remove exactly k jobs: P(removal = k)
                prob = dist.evalPMF(k);
            else
                % Remove all jobs (m = 1): P(removal >= n-1)
                cdfNMinus1 = 0;
                for j = 0:(n-2)
                    cdfNMinus1 = cdfNMinus1 + dist.evalPMF(j);
                end
                prob = 1 - cdfNMinus1;
            end
            if prob > 0
                L(n, m) = L(n, m) + rate * prob;
            end
        end
    end
end

% Single removal negative arrival transitions: n -> n-1 (only if queue non-empty)
if lambda_ir_neg_single > 0
    for n = 2:Np
        L(n, n-1) = L(n, n-1) + lambda_ir_neg_single;
    end
end

% Service rate at this station
mu_ir = sn.rates(ist, r);
if ~isnan(mu_ir) && mu_ir > 0
    % Departures to sink (use rtnodes with node indices)
    % Get the node index for this station
    nodeIdx = sn.stationToNode(ist);

    prob_sink = 0;
    if isfield(sn, 'rtnodes') && ~isempty(sn.rtnodes)
        nNodes = length(sn.nodetype);
        for jsnk = sinkNodes
            for s = 1:K
                % rtnodes indices: (nodeIdx-1)*K + classIdx
                fromIdx = (nodeIdx - 1) * K + r;
                toIdx = (jsnk - 1) * K + s;
                if fromIdx <= size(sn.rtnodes, 1) && toIdx <= size(sn.rtnodes, 2)
                    prob_sink = prob_sink + sn.rtnodes(fromIdx, toIdx);
                end
            end
        end
    end

    % Self-routing (stays at same station, same class)
    prob_self = rt((ist-1)*K + r, (ist-1)*K + r);

    % Combined local departure rate
    local_departure_rate = mu_ir * prob_sink;

    % Departure transitions: n -> n-1 (add to existing negative arrival effects)
    for n = 2:Np
        L(n, n-1) = L(n, n-1) + local_departure_rate;
    end

    % Self-service transitions: n -> n (diagonal, for phase transitions)
    for n = 2:Np
        L(n, n) = mu_ir * prob_self;
    end
end

end

function [QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N)
% RCAT_METRICS Convert RCAT solution to LINE performance metrics

M = sn.nstations;
K = sn.nclasses;

QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);

% Compute metrics for each (station, class) pair
for ist = 1:M
    for r = 1:K
        p = processMap(ist, r);
        if p > 0 && ~isempty(pi) && p <= length(pi) && ~isempty(pi{p})
            Np = length(pi{p});

            % Queue length: E[N] = sum_{n=0}^{Np-1} n * pi(n+1)
            QN(ist, r) = (0:(Np-1)) * pi{p}(:);

            % Utilization: P(N > 0) = 1 - pi(0) = 1 - pi{p}(1)
            UN(ist, r) = 1 - pi{p}(1);

            % Throughput: compute from utilization and service rate
            mu_ir = sn.rates(ist, r);
            if ~isnan(mu_ir) && mu_ir > 0
                TN(ist, r) = mu_ir * UN(ist, r);
            end
        end
    end
end

% Handle self-looping classes: they always stay at their reference station
% and share the server with other classes under PS scheduling.
% Only override if the method did not already compute SLC metrics (QN == 0).
if isfield(sn, 'isslc') && any(sn.isslc)
    for r = 1:K
        if sn.isslc(r)
            refst = sn.refstat(r);
            if refst > 0 && refst <= M && QN(refst, r) == 0
                % Self-looping class: all jobs stay at reference station
                QN(refst, r) = sn.njobs(r);

                % Service rate for this class
                mu_ir = sn.rates(refst, r);
                if ~isnan(mu_ir) && mu_ir > 0
                    nservers = sn.nservers(refst);
                    if isinf(nservers)
                        % Delay (infinite server): no capacity constraint,
                        % each job gets dedicated service
                        UN(refst, r) = QN(refst, r);
                        TN(refst, r) = mu_ir * QN(refst, r);
                    else
                        % Queue (finite server): capacity constraint applies
                        % Get utilization from other classes at this station
                        other_util = 0;
                        for s = 1:K
                            if s ~= r && ~sn.isslc(s)
                                other_util = other_util + UN(refst, s);
                            end
                        end

                        % Remaining capacity is shared with SLC
                        remaining_capacity = max(0, 1 - other_util);

                        % SLC utilization: min(demand, remaining capacity)
                        slc_demand = QN(refst, r) / mu_ir;
                        UN(refst, r) = min(slc_demand, remaining_capacity);
                        TN(refst, r) = mu_ir * UN(refst, r);
                    end
                end
            end
        end
    end
end

% Response times from Little's law: R = Q / T
for ist = 1:M
    for r = 1:K
        if TN(ist, r) > 0
            RN(ist, r) = QN(ist, r) / TN(ist, r);
        else
            RN(ist, r) = 0;
        end
    end
end

% System metrics
CN = zeros(1, K);
XN = zeros(1, K);

% For open classes: system throughput = arrival rate, system response time = sum of response times
for r = 1:K
    if sn.njobs(r) >= Inf  % Open class
        % System throughput equals arrival rate (from source)
        for ist = 1:M
            nodeIdx = sn.stationToNode(ist);
            if sn.nodetype(nodeIdx) == NodeType.Source
                XN(r) = sn.rates(ist, r);
                break;
            end
        end
        % System response time = sum over all stations
        CN(r) = sum(RN(:, r));
    else
        % Closed class: use reference station
        refst = sn.refstat(r);
        if refst > 0 && refst <= M
            XN(r) = TN(refst, r);
            if XN(r) > 0
                CN(r) = sn.njobs(r) / XN(r);
            end
        end
    end
end

end

function pi = birth_death_solve(Q)
% BIRTH_DEATH_SOLVE Solve equilibrium of a birth-death (tridiagonal) CTMC
%
% For a birth-death chain with birth rate lambda_n = Q(n, n+1) and
% death rate mu_n = Q(n, n-1), the equilibrium is computed using the
% recursion pi(n) = pi(n-1) * lambda(n-1) / mu(n).
%
% This is numerically stable and avoids the ill-conditioned linear system
% that plagues null-space methods for large state spaces.

n = size(Q, 1);
if n <= 1
    pi = 1;
    return;
end

pi = zeros(1, n);
pi(1) = 1.0;

for i = 2:n
    birth_rate = Q(i-1, i);
    death_rate = Q(i, i-1);
    if death_rate > 0
        pi(i) = pi(i-1) * birth_rate / death_rate;
    else
        pi(i) = 0;
    end
end

total = sum(pi);
if total > 0
    pi = pi / total;
else
    pi = ones(1, n) / n;
end

end

function result = is_tridiagonal(Q)
% IS_TRIDIAGONAL Check if a matrix is tridiagonal
n = size(Q, 1);
result = true;
for i = 1:n
    for j = 1:n
        if abs(i - j) > 1 && abs(Q(i, j)) > 1e-14
            result = false;
            return;
        end
    end
end
end
