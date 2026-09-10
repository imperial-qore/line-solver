function [QN,UN,RN,TN,CN,XN,runtime,method,tranSysState,tranSync,sn,QNCI,UNCI,RNCI,TNCI,ANCI,WNCI,StartN,PreemptN,tranStartTag,tranPreemptTag] = solver_ssa_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,...,STARTN,PREEMPTN,TRANSTARTTAG,TRANPREEMPTTAG] = SOLVER_SSA_ANALYZER(SN, OPTIONS)
%
% STARTN/PREEMPTN are the derived (station x class) service-start and
% preemption rates, and TRANSTARTTAG/TRANPREEMPTTAG the per-step tags of the
% sampled path. Both back-ends fill them, so the counters do not depend on
% which one the dispatch picks.
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

% A marking-dependent firing rate has no engine here; the predicate is the one
% SolverSSA.supportsModelMethod asks, so the report and this run agree.
[fdOk, fdWhy] = ssa_firingdep_refusal(sn);
if ~fdOk
    line_error(mfilename, fdWhy);
end

% Convert non-Markovian distributions to PH
% SSA draws a sample path and the fluid ODEs read mu*phi as a flow, so their
% surrogate must be a genuine phase-type: a matrix exponential has neither.
options.config.phfit = 'ph';
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
StartN = zeros(M,K);
PreemptN = zeros(M,K);
tranStartTag = {};
tranPreemptTag = {};

% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

LineConsole.loop('drawing the sample path: %g samples requested', options.samples);

% The per-feature NRM tests, read individually below because the explicit 'nrm'
% arm falls back to the serial engine with a warning that names the offending
% construct. SSA_NRM_ELIGIBLE and SSA_NRM_REFUSAL read the same body.
nrmGuards = ssa_nrm_guards(sn);

% -------------------------------------------------------------------------
% Pick analysis back-end
% -------------------------------------------------------------------------
switch options.method
    case {'default'}          % "default" prefers NRM for closed QNs with INF/PS
        % see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch)
        % A net the NRM's reaction builder cannot read (ssa_nrm_guards.spn)
        % falls through to the eligibility test below, which sends it to the
        % serial engine instead of letting the builder raise.
        if any(sn.nodetype == NodeType.Transition) && nrmGuards.gd && nrmGuards.spn
            line_debug('Default method: SPN detected, using NRM');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel);
            end
            runtime = toc(Tstart);
            return
        end
        % see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch)
        nrmSupported = ssa_nrm_eligible(sn);

        if nrmSupported
            line_debug('Default method: using NRM (Next Reaction Method)\n');
            line_debug('Using NRM method (fast path), calling solver_ssa_analyzer_nrm');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            % Compute CI using batch means if enabled
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel, options);
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
        if ~nrmGuards.renege
            % Phase-type patience would need the remaining-patience phase of
            % each waiting job, which the reaction network does not carry
            line_warning(mfilename, 'NRM supports only exponential (memoryless) patience for reneging; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~nrmGuards.balk
            % EXPECTED_WAIT / COMBINED balking depend on the mean waiting time,
            % which is not a function of the state vector
            line_warning(mfilename, 'NRM only supports QUEUE_LENGTH balking; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~nrmGuards.phase
            % Phase-type service is expanded exactly only at the INF/PS and
            % non-preemptive buffered families; a preemptive (LCFSPR) or polling
            % station with phase-type service needs the serial engine.
            line_warning(mfilename, 'NRM expands phase-type service only at INF/PS and non-preemptive buffered stations; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~nrmGuards.gd
            % A global (Whittle) dependence reads the whole population matrix,
            % which the per-station propensity closures do not receive
            line_warning(mfilename, 'NRM does not support a global dependence (setGlobalDependence); falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~nrmGuards.cache
            % A cache with a retrieval (delayed-hit) system sends a miss to a
            % queue and back, which the immediate class-switch cache model does
            % not yet reproduce; use the serial engine.
            line_warning(mfilename, 'NRM does not yet support the cache retrieval (delayed-hit) system; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~nrmGuards.spn
            % A non-exponential firing, a non-exponential or Place-less Source
            % arrival, or an infinite marking: the reaction builder has no form
            % for them, State.afterGlobalEvent does.
            line_warning(mfilename, 'NRM supports exponential firings and exponential Source arrivals into a Place only; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        else
            % Finite capacity regions run in the NRM under both rules (DROP
            % censors the refused transition, WAITQ parks it in a per-region
            % FIFO), so no region rule forces a serial fallback here.
            line_debug('Using explicit NRM method, calling solver_ssa_analyzer_nrm');

            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            sn.method = 'default/nrm';
            % Compute CI using batch means if enabled
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel, options);
            end
            runtime = toc(Tstart);
            return
        end
    case 'ssa'                      % alias for serial path below
        line_debug('Using ssa alias, redirecting to serial method');
        options.method = 'serial';
        sn.method = 'default/serial';
end

% SERIAL / PARALLEL ANALYSERS (legacy paths) ------------------------------
switch options.method
    case {'serial'}
        line_debug('Using serial method, calling solver_ssa_analyzer_serial');
        [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN,tranStartTag,tranPreemptTag] = ...
            solver_ssa_analyzer_serial(sn, init_state, options, false);
        method = 'serial';

    case {'para','parallel'}
        % Prefer the NRM on the same eligibility gate as the default path: an
        % NRM-eligible model runs on the fast single-run NRM rather than
        % replicated serial simulation.
        if ssa_nrm_eligible(sn)
            line_debug('Parallel method: model is NRM-eligible, using NRM');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel, options);
            end
            runtime = toc(Tstart);
            return
        end
        line_debug('Using parallel method, calling solver_ssa_analyzer_parallel');
        try
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN] = ...
                solver_ssa_analyzer_parallel(sn, init_state, options);
            method = 'parallel';
        catch ME
            if strcmp(ME.identifier,'MATLAB:spmd:NoPCT')
                line_printf(['Parallel Computing Toolbox unavailable – ',...
                    'falling back to serial SSA.\n']);
                [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn,StartN,PreemptN,tranStartTag,tranPreemptTag] = ...
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
    [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel, options);
end

runtime = toc(Tstart);
end

function [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel, options)
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

    % Discard initial transient before batch means: use
    % options.config.warmupfrac when set (> 0), else the legacy 10% discard
    if isfield(options,'config') && isfield(options.config,'warmupfrac') ...
            && ~isempty(options.config.warmupfrac) && options.config.warmupfrac > 0
        warmupfrac = options.config.warmupfrac;
    else
        warmupfrac = 0.1;
    end
    transientCutoff = max(1, floor(nSamples * warmupfrac));

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
