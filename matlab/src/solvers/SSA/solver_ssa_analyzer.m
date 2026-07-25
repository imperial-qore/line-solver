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

% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

% -------------------------------------------------------------------------
% Pick analysis back-end
% -------------------------------------------------------------------------
switch options.method
    case {'default'}          % "default" prefers NRM for closed QNs with INF/PS
        % see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch)
        if any(sn.nodetype == NodeType.Transition)
            line_debug('Default method: SPN detected, using NRM');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
                solver_ssa_analyzer_nrm(sn, options);
            method = 'nrm';
            if confintEnabled && ~isempty(tranSysState) && length(tranSysState) > 1
                [QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = ssa_compute_batch_means_ci(tranSysState, sn, confintLevel);
            end
            runtime = toc(Tstart);
            return
        end
        % see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch)
        nrmSupported = isNrmEligible(sn);

        if nrmSupported
            line_debug('Default method: using NRM (Next Reaction Method)\n');
            line_debug('Using NRM method (fast path), calling solver_ssa_analyzer_nrm');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
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
        if ~routingNrmOK(sn)
            % NRM routes departures via the static rt matrix (JSQ and memoryless
            % SQ are handled natively); the remaining state-dependent
            % routing strategies need the serial engine
            line_warning(mfilename, 'NRM does not support RL routing; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~renegeNrmOK(sn)
            % Phase-type patience would need the remaining-patience phase of
            % each waiting job, which the reaction network does not carry
            line_warning(mfilename, 'NRM supports only exponential (memoryless) patience for reneging; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~balkNrmOK(sn)
            % EXPECTED_WAIT / COMBINED balking depend on the mean waiting time,
            % which is not a function of the state vector
            line_warning(mfilename, 'NRM only supports QUEUE_LENGTH balking; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~phaseNrmOK(sn)
            % Phase-type service is expanded exactly only at the INF/PS and
            % non-preemptive buffered families; a preemptive (LCFSPR) or polling
            % station with phase-type service needs the serial engine.
            line_warning(mfilename, 'NRM expands phase-type service only at INF/PS and non-preemptive buffered stations; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        elseif ~cacheNrmOK(sn)
            % A cache with a retrieval (delayed-hit) system sends a miss to a
            % queue and back, which the immediate class-switch cache model does
            % not yet reproduce; use the serial engine.
            line_warning(mfilename, 'NRM does not yet support the cache retrieval (delayed-hit) system; falling back to the serial method.');
            options.method = 'serial';
            sn.method = 'default/serial';
        else
            % Finite capacity regions run in the NRM under both rules (DROP
            % censors the refused transition, WAITQ parks it in a per-region
            % FIFO), so no region rule forces a serial fallback here.
            line_debug('Using explicit NRM method, calling solver_ssa_analyzer_nrm');

            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
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
        [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
            solver_ssa_analyzer_serial(sn, init_state, options, false);
        method = 'serial';

    case {'para','parallel'}
        % Prefer the NRM on the same eligibility gate as the default path: an
        % NRM-eligible model runs on the fast single-run NRM rather than
        % replicated serial simulation.
        if isNrmEligible(sn)
            line_debug('Parallel method: model is NRM-eligible, using NRM');
            [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = ...
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

function ok = fcrNrmOK(sn)
% Finite capacity regions are supported under both rules. DROP destroys a
% refused job (the departure fires and the job never reaches the destination);
% WAITQ parks it in a per-region FIFO and admits it head-of-line as capacity
% frees. The NRM carries the FIFO explicitly (see fcrReleaseCascade), so no
% region rule forces a fallback. Linear-constraint and memory-budget regions
% ride the same admission test.
ok = true;
end

function ok = routingNrmOK(sn)
% True when every routing strategy in the model is one the NRM resolves at
% firing time. JSQ and SQ(d) select from the candidate queue lengths and
% RROBIN/WRROBIN walk a rotation pointer, which the NRM carries as auxiliary
% state alongside the buffers, since no rate reads it -- it only steers the
% destination draw, so it need not enter the reaction network. RL needs an
% external policy and still requires the serial engine.
ok = ~any(sn.routing(:) == RoutingStrategy.RL);
if ~ok
    return
end
end

function ok = balkNrmOK(sn)
% True when no station uses a balking strategy that the NRM cannot evaluate.
% QUEUE_LENGTH is a pure function of the state vector, so the NRM draws it at
% firing time; EXPECTED_WAIT and COMBINED depend on the mean waiting time and
% need the serial engine (State.afterEventStation rejects them likewise).
ok = true;
if ~isfield(sn,'balkingStrategy') || isempty(sn.balkingStrategy)
    return
end
bs = sn.balkingStrategy(:);
ok = all(bs == 0 | bs == BalkingStrategy.QUEUE_LENGTH);
end

function ok = cacheNrmOK(sn)
% True for every Cache node. The NRM models a cache access as an immediate
% state-dependent class switch (read -> hit/miss/retrieval) at the cache node,
% applying the same replacement logic as State.afterEventCache to the cache
% contents carried alongside the buffers, INCLUDING the retrieval (delayed-hit)
% system: a miss for an item not yet being fetched begins a retrieval (the job
% is routed to the fetch queue and returns to complete the miss), and a
% concurrent request for an item already being fetched is absorbed as a delayed
% hit -- matching the serial engine's sample-path semantics.
ok = true;
end

function ok = isNrmEligible(sn)
% True when the NRM engine can run this model. The scheduling list matches
% solver_ssa_analyzer_nrm's own validation (INF/PS family, non-preemptive
% buffered family, LCFSPR, PAS, POLLING); the per-feature guards exclude the
% sub-cases it cannot reproduce. The NRM simulates open and closed models alike
% -- it only lacks Fork/Join node handling -- so the gate is a Fork/Join
% exclusion, not the INF/PS-only sn_is_population_model. Used to prefer the NRM
% on the default and parallel dispatch paths.
allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
    SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS, ...
    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
    SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING];
ok = all(arrayfun(@(s) any(s == allowedSched), sn.sched)) && ...
    cacheNrmOK(sn) && ...
    ~sn_has_fork_join(sn) && ...
    routingNrmOK(sn) && ...
    fcrNrmOK(sn) && ...
    balkNrmOK(sn) && ...
    renegeNrmOK(sn) && ...
    phaseNrmOK(sn);
end

function ok = renegeNrmOK(sn)
% True when no station renegs with non-exponential patience. The NRM abandons
% at the aggregate rate (waiting count)*mu, which is only correct when patience
% is memoryless; phase-type patience would need each waiting job's remaining
% phase. SOLVER_SSA rejects the same combination outright.
ok = true;
if ~isfield(sn,'impatienceClass') || isempty(sn.impatienceClass)
    return
end
bad = (sn.impatienceClass == ImpatienceType.RENEGING) & (sn.impatienceType ~= ProcessType.EXP);
ok = ~any(bad(:));
end

function ok = phaseNrmOK(sn)
% True when every non-exponential service sits at a station whose rate law the
% NRM expands exactly. Phase expansion splits the class-level share across a
% class's phases in the ratio kir/nir, which needs only the per-phase
% populations -- true of the INF/PS family, where every job present is in
% service. A buffered policy instead needs the phase multiset of the jobs
% ACTUALLY in service, which the waiting-only buffer does not record, so
% non-exponential service there still needs the serial engine.
ok = true;
% see _kb/06-solver-catalog.md for rationale (SSA NRM dispatch, EXT-source exclusion)
exact = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.LPS, ...
    SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, ...
    SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT];
for ist = 1:sn.nstations
    for r = 1:sn.nclasses
        if sn.procid(ist,r) == ProcessType.DISABLED || sn.procid(ist,r) == ProcessType.EXP
            continue
        end
        if ~any(sn.sched(ist) == exact)
            ok = false;
            return
        end
    end
end
end
