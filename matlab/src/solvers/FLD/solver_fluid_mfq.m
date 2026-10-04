function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_mfq(sn, options)
%SOLVER_FLUID_MFQ Solves single-queue open systems using BUTools MAPMAP1
%
% [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = SOLVER_FLUID_MFQ(sn, options)
%
% Uses the BUTools library's FluFluQueue function to analyze MAP/MAP/1 queues.
% Provides exact steady-state queue length and sojourn time moments for
% single-queue open systems with phase-type arrivals and service.
%
% This solver applies only to single-queue topologies:
%   Source (external arrivals) -> Queue (finite server) -> Sink
%
% Applicability:
%   - Single-queue open system
%   - Single-server (c=1) or infinite-server (c=Inf)
%   - No feedback loops
%   - Supported distributions: Exp, Erlang, HyperExp, Cox, APH, MAP, PH, MMDP
%
% Parameters:
%   sn (struct): Network structure from Network.getStruct()
%   options (struct): Solver options including:
%       - method: 'mfq' (should be set by dispatcher)
%       - tol: Numerical tolerance (default: 1e-14)
%
% Returns:
%   QN (MxK): Mean queue lengths at each station
%   UN (MxK): Utilizations at each station
%   RN (MxK): Mean response times at each station
%   TN (MxK): Throughputs at each station
%   xvec_it (cell): Final state vector for compatibility
%   QNt (MxK cell): Transient queue lengths (empty for MAPMAP1)
%   UNt (MxK cell): Transient utilizations (empty for MAPMAP1)
%   TNt (MxK cell): Transient throughputs (empty for MAPMAP1)
%   xvec_t: Transient state vectors (empty for MAPMAP1)
%   t: Time vector
%   iters: Number of iterations (always 1 for MAPMAP1)
%   runtime: Execution time in seconds
%
% References:
%   BUTools: Queueing and traffic modeling library
%   https://github.com/ghorvath78/butools
%
% See also: solver_fluid_analyzer, fluid_is_single_queue, FluFluQueue

    runtime_start = tic;

    M = sn.nstations;
    K = sn.nclasses;

    % Initialize output arrays
    QN = zeros(M, K);
    UN = zeros(M, K);
    RN = zeros(M, K);
    TN = zeros(M, K);

    % Initialize transient outputs (MFQ is steady-state only)
    t = [0; options.timespan(2)];
    QNt = cell(M, K);
    UNt = cell(M, K);
    TNt = cell(M, K);
    for ist = 1:M
        for k = 1:K
            QNt{ist, k} = zeros(2, 1);
            UNt{ist, k} = zeros(2, 1);
            TNt{ist, k} = zeros(2, 1);
        end
    end

    xvec_t = [];
    xvec_it = {zeros(size(sn.state{1}, 2), 1)};

    % =========================================================================
    % TOPOLOGY VALIDATION
    % =========================================================================

    [isSingleQueue, fluidInfo] = fluid_is_single_queue(sn);

    if ~isSingleQueue
        line_error(mfilename, 'MFQ requires single-queue topology: %s', fluidInfo.errorMsg);
    end

    % Get station indices (for accessing sn arrays)
    sourceIdx = fluidInfo.sourceStation;
    queueIdx = fluidInfo.queueStation;
    % Note: sinkIdx is a node index, but we don't need to access station arrays for sink

    line_debug('MFQ: Source station=%d, Queue station=%d', sourceIdx, queueIdx);

    % =========================================================================
    % EXTRACT PARAMETERS FOR EACH CLASS
    % =========================================================================

    % Process each class independently
    for k = 1:K
        if isinf(sn.njobs(k))
            % Open class
            line_debug('MFQ: Processing open class %d', k);

            % Extract MFQ parameters
            [Qin, Rin, Qout, Rout, srv0stop, fluidInfo_k] = fluid_extract_params(sn, fluidInfo);

            % Get arrival rate
            lambda = sn.rates(sourceIdx, k);

            % ===================================================================
            % DETERMINE IF SIMPLE EXPONENTIAL (M/M/1) OR COMPLEX (MAP/PH)
            % ===================================================================

            % Check if both arrival and service are simple exponential (1-state)
            is_simple_exponential = isscalar(Qin) && isscalar(Qout) && Qin == 0 && Qout == 0;

            if is_simple_exponential
                % Simple M/M/1 case: use analytical formulas
                % FluFluQueue doesn't handle trivial 1-state CTMCs well
                line_debug('MFQ: Using analytical M/M/1 formulas for class %d', k);

                mu = fluidInfo_k.mu;
                rho = lambda / mu;

                if rho >= 1
                    line_warning(mfilename, 'System is unstable (rho=%.4f >= 1) for class %d', rho, k);
                    flMoms = [Inf, Inf];
                    stMoms = [Inf, Inf];
                else
                    % M/M/1 analytical formulas
                    % E[L] = rho / (1 - rho)
                    % E[W] = 1 / (mu - lambda)
                    % Var[L] = rho / (1 - rho)^2
                    % Var[W] = 1 / (mu - lambda)^2
                    L_mean = rho / (1 - rho);
                    L_var = rho / (1 - rho)^2;
                    W_mean = 1 / (mu - lambda);
                    W_var = 1 / (mu - lambda)^2;

                    flMoms = [L_mean, L_mean^2 + L_var];  % E[L], E[L^2]
                    stMoms = [W_mean, W_mean^2 + W_var];  % E[W], E[W^2]
                end
            else
                % Complex case: call FluFluQueue
                try
                    % Request 2 moments for mean and variance
                    % Set precision to match solver tolerance
                    prec = max(options.tol, 1e-14);

                    line_debug('MFQ: Calling FluFluQueue for class %d with precision %.2e', k, prec);

                    % Call FluFluQueue to get fluid level moments and sojourn time moments
                    [flMoms, stMoms] = FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, ...
                        'flMoms', 2, 'stMoms', 2, 'prec', prec);
                catch ME
                    line_error(mfilename, 'FluFluQueue failed for class %d: %s', k, ME.message);
                end
            end

            % ===============================================================
            % MAP RESULTS TO LINE METRICS
            % ===============================================================

            % Fluid level first moment = mean queue length
            QN(queueIdx, k) = flMoms(1);

            % Sojourn time first moment = mean response time
            RN(queueIdx, k) = stMoms(1);

            % Throughput from Little's Law: X = L / W
            if stMoms(1) > GlobalConstants.FineTol
                TN(queueIdx, k) = flMoms(1) / stMoms(1);
            else
                TN(queueIdx, k) = lambda;
            end

            % Utilization: U = lambda * E[S]
            % Extract mean service time from service process
            mu = fluidInfo_k.mu;
            UN(queueIdx, k) = min(1, lambda * (1 / mu));

            % Set throughput at source equal to queue throughput
            TN(sourceIdx, k) = TN(queueIdx, k);

            % ===============================================================
            % TRANSIENT SOLUTION (TRIVIAL FOR MFQ)
            % ===============================================================

            % MFQ provides steady-state solution only
            % Return trivial transient: jump from 0 to steady-state

            QNt{queueIdx, k} = [0; QN(queueIdx, k)];
            UNt{queueIdx, k} = [0; UN(queueIdx, k)];
            TNt{queueIdx, k} = [0; TN(queueIdx, k)];

            QNt{sourceIdx, k} = [0; QN(sourceIdx, k)];
            UNt{sourceIdx, k} = [0; UN(sourceIdx, k)];
            TNt{sourceIdx, k} = [0; TN(sourceIdx, k)];

            % ===============================================================
            % VALIDATION: LITTLE'S LAW
            % ===============================================================

            little_lhs = flMoms(1);
            little_rhs = lambda * stMoms(1);

            if little_rhs > 0
                little_error = abs(little_lhs - little_rhs) / little_rhs;
                if little_error > 0.01  % 1% tolerance
                    line_warning(mfilename, ...
                        'Little''s Law violation for class %d: L=%.6f vs λW=%.6f (error=%.2f%%)', ...
                        k, little_lhs, little_rhs, little_error * 100);
                else
                    line_debug('Little''s Law satisfied for class %d (error=%.4f%%)', k, little_error * 100);
                end
            end

        else
            % Closed class: not supported by MFQ
            line_warning(mfilename, 'MFQ only supports open classes. Class %d is closed, skipping.', k);
        end
    end

    % =========================================================================
    % FINALIZE OUTPUT
    % =========================================================================

    % System-level metrics (from reference station)
    XN = zeros(1, K);
    CN = zeros(1, K);

    for k = 1:K
        if sn.refstat(k) > 0  % Ignore artificial classes
            XN(k) = TN(sn.refstat(k), k);
            if isinf(sn.njobs(k))
                % Open class: CN = sum of response times
                CN(k) = sum(RN(:, k));
            else
                % Closed class: CN = N / X (should not reach here)
                CN(k) = sn.njobs(k) / XN(k);
            end
        end
    end

    % Assign system metrics to source station for compatibility
    % (source is the reference station for open classes)
    for k = 1:K
        if isinf(sn.njobs(k))
            TN(sourceIdx, k) = XN(k);
        end
    end

    % Create state vector (dummy for compatibility)
    xvec_it = {zeros(1, 1)};

    % Number of iterations is 1 for MFQ (no iteration needed)
    iters = 1;

    runtime = toc(runtime_start);

    line_debug('MFQ completed in %.4f seconds', runtime);
end

function [Qin, Rin, Qout, Rout, srv0stop, fluidInfo] = fluid_extract_params(sn, fluidInfo)
%FLUID_EXTRACT_PARAMS Extract FluFluQueue parameters from LINE network structure

    % Use station indices (not node indices)
    sourceIdx = fluidInfo.sourceStation;
    queueIdx = fluidInfo.queueStation;

    % For now, analyze single-class (class 1)
    classIdx = 1;

    % Default boundary behavior: service stops when empty (work-conserving)
    srv0stop = true;

    % Override if specified in options
    if isfield(fluidInfo, 'srv0stop')
        srv0stop = fluidInfo.srv0stop;
    end

    % =========================================================================
    % EXTRACT ARRIVAL PROCESS (from Source station)
    % =========================================================================

    % Get arrival rate
    lambda = sn.rates(sourceIdx, classIdx);

    if isnan(lambda) || lambda <= 0
        error('No valid arrival rate for class %d at source', classIdx);
    end

    % Get arrival process
    arrivalProc = sn.proc{sourceIdx}{classIdx};

    if isempty(arrivalProc)
        % Simple Poisson arrival (single state)
        Qin = 0;
        Rin = lambda;
    elseif isa(arrivalProc, 'MMDP')
        % MMDP: use Q and R directly (already in BUTools format)
        Qin = arrivalProc.Q();
        Rin = arrivalProc.R();
    else
        % MAP or PH arrival process
        % Check if it's a simple exponential (1 state)
        if size(arrivalProc{1}, 1) == 1
            % Single state: Poisson arrival
            Qin = 0;
            Rin = lambda;
        else
            % Multi-state MAP/MMPP: convert to CTMC format
            [Qin, Rin] = fluid_dist2butools(arrivalProc, lambda, true);
        end
    end

    % =========================================================================
    % EXTRACT SERVICE PROCESS (from Queue station)
    % =========================================================================

    % Get service rates per phase
    mu_vec = sn.mu{queueIdx}{classIdx};

    if isempty(mu_vec)
        error('No valid service rate for class %d at queue', classIdx);
    end

    % Average service rate (for normalization)
    mu_avg = mean(mu_vec);

    % Get service process
    serviceProc = sn.proc{queueIdx}{classIdx};

    if isempty(serviceProc)
        % Simple exponential service (single phase)
        Qout = 0;
        Rout = mu_avg;
    elseif isa(serviceProc, 'MMDP')
        % MMDP: use Q and R directly (already in BUTools format)
        Qout = serviceProc.Q();
        Rout = serviceProc.R();
    else
        % PH or Coxian service process
        % Check if it's a simple exponential (1 phase)
        if size(serviceProc{1}, 1) == 1
            % Single phase: exponential service
            Qout = 0;
            Rout = mu_avg;
        else
            % Multi-phase: general phase-type service
            [Qout, Rout] = fluid_dist2butools(serviceProc, mu_avg, false);
        end
    end

    % =========================================================================
    % VALIDATION
    % =========================================================================

    % Validate dimensions
    if isscalar(Qin)
        % Trivial case: scalar CTMC (Poisson arrival)
        if Qin ~= 0 || Rin <= 0
            error('Invalid arrival process parameters');
        end
    else
        [n_in, m_in] = size(Qin);
        if n_in ~= m_in
            error('Qin must be square');
        end
        [n_rin, m_rin] = size(Rin);
        if n_rin ~= n_in || m_rin ~= n_in
            error('Rin dimensions do not match Qin');
        end
    end

    if isscalar(Qout)
        % Trivial case: scalar CTMC (exponential service)
        if Qout ~= 0 || Rout <= 0
            error('Invalid service process parameters');
        end
    else
        [n_out, m_out] = size(Qout);
        if n_out ~= m_out
            error('Qout must be square');
        end
        [n_rout, m_rout] = size(Rout);
        if n_rout ~= n_out || m_rout ~= n_out
            error('Rout dimensions do not match Qout');
        end
    end

    % Store extracted parameters in fluidInfo for later use
    fluidInfo.lambda = lambda;
    fluidInfo.mu = mu_avg;
    fluidInfo.arrivalProc = arrivalProc;
    fluidInfo.serviceProc = serviceProc;
end

function [Q, R] = fluid_dist2butools(ph_process, rate_scale, is_arrival)
%FLUID_DIST2BUTOOLS Convert LINE phase-type distribution to BUTools CTMC format

    % Extract phase-type matrices
    if ~iscell(ph_process) || length(ph_process) < 2
        error('ph_process must be a cell array {D0, D1}');
    end

    D0 = ph_process{1};
    D1 = ph_process{2};

    % Check dimensions
    [n_d0, m_d0] = size(D0);
    [n_d1, m_d1] = size(D1);

    if n_d0 ~= m_d0
        error('D0 must be square');
    end

    if n_d1 ~= n_d0 || m_d1 ~= m_d0
        error('D0 and D1 must have the same dimensions');
    end

    n = n_d0;  % Number of states

    % Compute generator Q = D0 + D1
    Q = D0 + D1;

    % Validate Q is a proper generator
    for i = 1:n
        % Diagonal element should be negative or zero
        if Q(i, i) > 1e-14
            warning('Generator Q has positive diagonal element Q(%d,%d)=%.4e', i, i, Q(i, i));
        end

        % Off-diagonal elements should be non-negative
        for j = 1:n
            if i ~= j && Q(i, j) < -1e-14
                warning('Generator Q has negative off-diagonal element Q(%d,%d)=%.4e', i, j, Q(i, j));
            end
        end

        % Row sum should be non-positive (≤ 0)
        rowsum = sum(Q(i, :));
        if rowsum > 1e-14
            warning('Generator Q has positive row sum at row %d: sum=%.4e', i, rowsum);
        end
    end

    % Compute rate matrix R
    R = zeros(n, n);

    % For both arrival and service, use D1 to determine fluid rates
    for i = 1:n
        R(i, i) = sum(D1(i, :));
    end

    % Validate rate matrix
    if any(diag(R) < -1e-14)
        error('Rate matrix R has negative diagonal entries');
    end

    % Check if R is actually diagonal (off-diagonal should be zero)
    off_diag_sum = sum(sum(abs(R - diag(diag(R)))));
    if off_diag_sum > 1e-12
        warning('Rate matrix R is not diagonal, taking diagonal part only');
        R = diag(diag(R));
    end
end
