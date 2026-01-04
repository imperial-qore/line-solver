function [QN,UN,RN,TN,CN,XN,totiter,percResults] = solver_mam_fj(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER,PERCRESULTS] = SOLVER_MAM_FJ(SN, OPTIONS)
%
% Analyze Fork-Join network using FJ_codes percentile approximation
%
% Parameters:
%   sn - Network structure from getStruct()
%   options - Solver options with FJ-specific fields:
%             - percentiles: percentile levels (default [0.90, 0.95, 0.99])
%             - fj_accuracy: C parameter for approximation (default 100)
%             - fj_tmode: 'NARE' or 'Sylves' for T matrix computation (default 'NARE')
%
% Returns:
%   QN, UN, RN, TN, CN, XN - Standard LINE performance metrics matrices
%   totiter - Number of iterations (N/A for FJ, returns 0)
%   percResults - Percentile results structure with fields:
%                 .RT - cell array of percentile structs per class
%                 .K - number of parallel queues
%                 .method - 'qiu'
%
% Reference:
%   Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
%   Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
%   Copyright 2015 Imperial College London
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate FJ topology
[isFJ, fjInfo] = fj_isfj(sn);

if ~isFJ
    line_error(mfilename, 'Network is not a valid Fork-Join topology for FJ_codes: %s', fjInfo.errorMsg);
end

% Extract FJ parameters
[arrival, service, K, fjInfo] = fj_extract_params(sn, fjInfo);

% Get FJ-specific options from config
% Request dense percentiles for accurate mean calculation via numerical integration
% E[X] = integral_0^inf (1 - F(x)) dx, approximated using trapezoidal rule
pers_dense = [0.01:0.05:0.95, 0.99, 0.999]; % Dense percentiles for mean calculation
pers_stored = [0.50, 0.90, 0.95, 0.99]; % Percentiles stored for user retrieval

if isfield(options.config, 'fj_accuracy')
    Cs = options.config.fj_accuracy;
else
    Cs = 100; % Default accuracy parameter
end

if isfield(options.config, 'fj_tmode')
    T_Mode = options.config.fj_tmode;
else
    T_Mode = 'NARE'; % Default method for T matrix
end

% Initialize result matrices
M = sn.nstations;
R = sn.nclasses;
QN = zeros(M, R);
UN = zeros(M, R);
RN = zeros(M, R);
TN = zeros(M, R);
CN = zeros(M, R);
XN = zeros(M, R);

% Initialize percentile results
percResults = struct();
percResults.RT = cell(1, R);
percResults.K = K;
percResults.method = 'qiu';
percResults.percentiles = pers_stored;
percResults.accuracy = Cs;
percResults.tmode = T_Mode;

% Analyze each class independently
for r = 1:R
    % Call FJ_codes mainFJ function with dense percentiles for mean calculation
    percentileRT_dense = mainFJ(arrival{r}, service{r}, pers_dense, K, Cs, T_Mode);

    % Call FJ_codes again with stored percentiles for user retrieval
    percentileRT_stored = mainFJ(arrival{r}, service{r}, pers_stored, K, Cs, T_Mode);

    % Store percentile results for this class (stored percentiles only)
    percResults.RT{r} = percentileRT_stored{1}; % mainFJ returns cell array, get first element

    % ===== Compute mean response time from percentile distribution =====
    % Use numerical integration: E[X] = integral_0^inf (1 - F(x)) dx
    % Approximated using trapezoidal rule on the inverse CDF (percentile function)

    % Extract percentile values and probabilities
    percentile_probs = percentileRT_dense{1}.percentiles / 100; % Convert to [0,1]
    percentile_values = percentileRT_dense{1}.RTp;

    % Add boundary points for better integration
    % At p=0, x=0 (assuming response time starts at 0)
    % At p=1, extrapolate or use last value
    probs_extended = [0; percentile_probs(:); 1];
    values_extended = [0; percentile_values(:); percentile_values(end) * 1.1]; % Small extrapolation at tail

    % Compute mean using trapezoidal integration of inverse CDF
    % E[X] = integral_0^1 Q(p) dp, where Q(p) is the inverse CDF (percentile function)
    mean_fj_rt = trapz(probs_extended, values_extended);

    % ===== Compute average metrics from FJ_codes results =====

    % Get arrival rate and service rate
    lambda = arrival{r}.lambda;
    mu = service{r}.mu;

    % For Fork-Join, the effective service time is the max of K parallel servers
    % We can use the percentile data to estimate mean response time
    % Or compute it analytically for exponential case

    % Mean response time for FJ queue (approximation using first-order moment)
    % For exponential service: E[R] ≈ sum(1/(mu*k)) for k=1 to K
    if service{r}.SerChoice == 1 % Exponential
        mean_service_time = sum(1./(mu * (1:K)));
    else
        % For general distributions, use the 50th percentile as rough estimate
        % or compute from distribution moments
        mean_service_time = 1/mu * sum(1./(1:K)); % Rough approximation
    end

    % Get station indices for queues between fork and join
    queueIdx = fjInfo.queueIdx;
    forkIdx = fjInfo.forkIdx;
    joinIdx = fjInfo.joinIdx;
    sourceIdx = find(sn.nodetype == NodeType.Source, 1);
    sinkIdx = find(sn.nodetype == NodeType.Sink, 1);

    % Throughput: arrival rate (assuming stable system)
    throughput = lambda;

    % Compute metrics for each queue
    for k = 1:K
        queueStat = sn.nodeToStation(queueIdx(k));

        % Each queue sees arrival rate lambda and has service rate mu
        % Utilization: rho = lambda/mu
        rho_k = lambda / mu;
        UN(queueStat, r) = rho_k;

        % Throughput at each queue
        TN(queueStat, r) = lambda;

        % Mean queue length: Using M/M/1 approximation
        QN(queueStat, r) = rho_k / (1 - rho_k);

        % Mean response time at queue: Using M/M/1 formula
        RN(queueStat, r) = 1 / (mu - lambda);
    end

    % Fork node (no queueing)
    if ~isnan(sn.nodeToStation(forkIdx))
        forkStat = sn.nodeToStation(forkIdx);
        TN(forkStat, r) = throughput;
        UN(forkStat, r) = 0;
        QN(forkStat, r) = 0;
        RN(forkStat, r) = 0;
    end

    % Join node (synchronization delay captured in FJ analysis)
    if ~isnan(sn.nodeToStation(joinIdx))
        joinStat = sn.nodeToStation(joinIdx);
        TN(joinStat, r) = throughput;
        UN(joinStat, r) = 0;

        % FJ_codes computes the TOTAL Fork-Join response time (from fork to join completion).
        % This includes both the individual queue service time AND the synchronization delay.
        % To report Join's contribution separately, we subtract the mean individual queue RT.
        mean_individual_queue_rt = 1 / (mu - lambda);
        synchronization_delay = mean_fj_rt - mean_individual_queue_rt;

        % Mean response time at Join: synchronization overhead only
        % This represents the additional delay due to waiting for the slowest of K parallel tasks
        RN(joinStat, r) = synchronization_delay;

        % Mean queue length at Join using Little's Law: L = lambda * W
        QN(joinStat, r) = throughput * synchronization_delay;
    end

    % Source node
    if ~isnan(sn.nodeToStation(sourceIdx))
        sourceStat = sn.nodeToStation(sourceIdx);
        TN(sourceStat, r) = throughput;
        UN(sourceStat, r) = 0;
        QN(sourceStat, r) = 0;
        RN(sourceStat, r) = 0;
    end

    % Sink node
    if ~isnan(sn.nodeToStation(sinkIdx))
        sinkStat = sn.nodeToStation(sinkIdx);
        TN(sinkStat, r) = throughput;
        UN(sinkStat, r) = 0;
        QN(sinkStat, r) = 0;
        RN(sinkStat, r) = 0;
    end
end

% No iterations for FJ_codes (it's an analytical/approximation method)
totiter = 0;

end
