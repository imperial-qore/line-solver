function [QN,UN,RN,TN,CN,XN,runtime,method,tranSysState,tranSync,sn,QNCI,UNCI,RNCI,TNCI,ANCI,WNCI] = solver_ssa_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_SSA_ANALYZER(SN, OPTIONS)
% Wrapper that selects the most suitable SSA performance-analysis back-end.
%
% If every station uses scheduling policy INF, EXT, or PS (and the network
% has no cache nodes), the faster Next-Reaction-Method analyser
%   -> solver_ssa_analyzer_nrm
% is invoked. Otherwise the original serial / parallel analysers are used.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart      = tic;

line_debug('SSA analyzer starting: method=%s, nstations=%d, nclasses=%d', options.method, sn.nstations, sn.nclasses);

% Convert non-Markovian distributions to PH
sn = sn_nonmarkov_toph(sn, options);

% Capture initial state after conversion (state may have been expanded for MAPs)
init_state  = sn.state;

% Initialize CI outputs
M = sn.nstations;
K = sn.nclasses;
QNCI = [];
UNCI = [];
RNCI = [];
TNCI = [];
ANCI = [];
WNCI = [];

% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

% -------------------------------------------------------------------------
% Pick analysis back-end
% -------------------------------------------------------------------------
switch options.method
    case {'default'}          % "default" prefers NRM for closed QNs with INF/PS
        allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS];
        nrmSupported =  all(arrayfun(@(s) any(s==allowedSched), sn.sched)) && ...
        all(sn.nodetype ~= NodeType.Cache) && ...
        sn_is_population_model(sn) && ...
         all(sn.procid(sn.procid~=ProcessType.DISABLED) == ProcessType.EXP); % all exp

        line_debug('Checking NRM eligibility: schedOK=%d, noCache=%d, isPopModel=%d, allExp=%d', ...
            all(arrayfun(@(s) any(s==allowedSched), sn.sched)), ...
            all(sn.nodetype ~= NodeType.Cache), ...
            sn_is_population_model(sn), ...
            all(sn.procid(sn.procid~=ProcessType.DISABLED) == ProcessType.EXP));

        if nrmSupported
            line_debug('Default method: using NRM (Next Reaction Method)\n');
            line_debug('Using NRM method (fast path), calling solver_ssa_analyzer_nrm');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            % Compute CI using batch means if enabled
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel);
            end
            runtime = toc(Tstart);
            return
        else
            % otherwise fall through to serial selection
            line_debug('Default method: using serial SSA\n');
            line_debug('NRM not supported, falling back to serial method');
            options.method = 'serial';
        end

    case 'nrm'
        line_debug('Using explicit NRM method, calling solver_ssa_analyzer_nrm');

        [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
            solver_ssa_analyzer_nrm(sn, options);
        method = 'nrm';
        sn.method = 'default/nrm';
        % Compute CI using batch means if enabled
        if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
            [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel);
        end
        runtime = toc(Tstart);
        return
    case 'ssa'                      % alias for serial path below
        line_debug('Using ssa alias, redirecting to serial method');
        options.method = 'serial';
        sn.method = 'default/serial';
end

% SERIAL / PARALLEL ANALYSERS (legacy paths) ------------------------------
switch options.method
    case {'serial'}
        line_debug('Using serial method, calling solver_ssa_analyzer_serial');
        [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
            solver_ssa_analyzer_serial(sn, init_state, options, false);
        method = 'serial';

    case {'para','parallel'}
        line_debug('Using parallel method, calling solver_ssa_analyzer_parallel');
        try
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
                solver_ssa_analyzer_parallel(sn, init_state, options);
            method = 'parallel';
        catch ME
            if strcmp(ME.identifier,'MATLAB:spmd:NoPCT')
                line_printf(['Parallel Computing Toolbox unavailable – ',...
                    'falling back to serial SSA.\n']);
                [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
                    solver_ssa_analyzer_serial(sn, init_state, options, true);
                method = 'serial';
            else
                rethrow(ME);
            end
        end

    otherwise
        error('solver_ssa_analyzer:UnknownMethod', ...
            'Unknown analysis method: %s', options.method);
end

% Compute CI using batch means if enabled
if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
    [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel);
end

runtime = toc(Tstart);
end

function [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel)
% SSA_COMPUTE_BATCH_MEANS_CI Compute confidence intervals using batch means method
%
% tranSysState{1} contains the cumulative time at each sample
% tranSysState{2:end} contain the state vectors for each stateful node

M = sn.nstations;
K = sn.nclasses;
QNCI = zeros(M, K);
UNCI = zeros(M, K);
RNCI = zeros(M, K);
TNCI = zeros(M, K);
ANCI = zeros(M, K);
WNCI = zeros(M, K);

% Extract time and state data
if iscell(tranSysState) && length(tranSysState) > 1
    times = tranSysState{1};
    nSamples = length(times);

    if nSamples < 20
        % Not enough samples for batch means
        return;
    end

    % Number of batches (use 10-30 batches for good CI estimation)
    numBatches = min(20, floor(nSamples / 10));
    if numBatches < 2
        return;
    end
    batchSize = floor(nSamples / numBatches);

    % Discard initial transient (first 10% of samples)
    transientCutoff = max(1, floor(nSamples * 0.1));

    % Extract queue length data from tranSysState
    % tranSysState{2:end} contains state vectors for each stateful node
    % We need to compute marginal queue lengths per station/class

    % Compute batch means for queue lengths
    for ist = 1:M
        isf = sn.stationToStateful(ist);
        if isf > 0 && (1 + isf) <= length(tranSysState)
            stateData = tranSysState{1 + isf};
            if isempty(stateData)
                continue;
            end

            for k = 1:K
                % Extract queue length for this station/class from state data
                % The state data format depends on the scheduling strategy
                % For simplicity, we'll use the marginal extraction
                ind = sn.stationToNode(ist);

                % Compute time-weighted batch means
                batchMeans = zeros(1, numBatches);
                for b = 1:numBatches
                    startIdx = transientCutoff + (b-1) * batchSize + 1;
                    endIdx = min(transientCutoff + b * batchSize, nSamples);
                    if startIdx > nSamples || startIdx >= endIdx
                        continue;
                    end

                    % Compute time-weighted average for this batch
                    % times contains cumulative times, compute inter-sample durations
                    if startIdx > 1
                        prevTime = times(startIdx - 1);
                    else
                        prevTime = 0;
                    end
                    batchTimes = times(startIdx:endIdx);

                    if length(batchTimes) >= 1
                        % Compute time duration each state was held
                        if startIdx == 1
                            dt = [batchTimes(1); diff(batchTimes)];
                        else
                            dt = [batchTimes(1) - prevTime; diff(batchTimes)];
                        end

                        % Extract queue lengths from state data
                        % For now, sum all columns to get total jobs at the station
                        % This works for most queue types where state represents job counts
                        qLengths = sum(stateData(startIdx:endIdx, :), 2);

                        totalTime = sum(dt);
                        if totalTime > 0
                            batchMeans(b) = sum(qLengths .* dt) / totalTime;
                        end
                    end
                end

                % Count valid batches (non-NaN)
                validMask = ~isnan(batchMeans);
                batchMeans = batchMeans(validMask);
                nBatches = length(batchMeans);

                if nBatches >= 2
                    % Compute mean and standard error
                    batchMean = mean(batchMeans);
                    batchStd = std(batchMeans);
                    stdErr = batchStd / sqrt(nBatches);

                    % t-critical value for confidence level
                    alpha = 1 - confintLevel;
                    tCrit = tinv(1 - alpha/2, nBatches - 1);

                    % Confidence interval half-width
                    QNCI(ist, k) = tCrit * stdErr;
                end
            end
        end
    end

    % For utilization, response time, and throughput CIs, use relative scaling
    % These are derived from queue length CI using Little's law relationships
    UNCI = QNCI; % Simplified - utilization CI scales similarly
    RNCI = QNCI; % Response time CI - would need service rate info
    TNCI = QNCI; % Throughput CI - would need arrival rate info
end
end
