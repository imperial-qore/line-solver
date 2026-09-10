function [Qt, Ut, Tt] = solver_mam_ldqbd_transient(sn, options)
% SOLVER_MAM_LDQBD_TRANSIENT Transient analysis of open queues via standard QBD
%
% Supports single-class open models (Source -> Queue -> Sink).
% For M/M/c: scalar QBD levels (any number of servers).
% For M/PH/1: m-phase QBD levels (single server only).
%
% Infinite capacity: libQBD adaptive Taylor series.
% Finite capacity: direct matrix exponentiation (expm).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Validate model structure
M = sn.nstations;
K = sn.nclasses;

if K ~= 1
    line_error(mfilename, 'Transient QBD method requires a single-class model.');
end

N = sn.njobs';
if ~isinf(N)
    line_error(mfilename, 'Transient QBD method requires an open model.');
end

%% Identify stations
sourceIdx = find(sn.sched == SchedStrategy.EXT);
queueIdx = find(sn.sched == SchedStrategy.FCFS);

if numel(sourceIdx) ~= 1 || numel(queueIdx) ~= 1
    line_error(mfilename, 'Transient QBD method requires exactly one Source and one FCFS Queue.');
end

%% Extract parameters
lambda = sn.rates(sourceIdx, 1);
PH_queue = sn.proc{queueIdx}{1};
nServers = sn.nservers(queueIdx);
bufCap = sn.cap(queueIdx);

if numel(PH_queue{1}) == 1
    mu = -PH_queue{1};
    nPhases = 1;
    isPH = false;
else
    nPhases = size(PH_queue{1}, 1);
    isPH = true;
    D0 = PH_queue{1};
    D1 = PH_queue{2};
    alpha = map_pie(PH_queue);
    t_exit = -D0 * ones(nPhases, 1);
end

if isPH && nServers > 1
    line_error(mfilename, 'Transient QBD with PH service supports single-server only.');
end

%% Time parameters
T_start = options.timespan(1);
T_end = options.timespan(2);
T_duration = T_end - T_start;

if ~isinf(bufCap)
    %% Finite capacity: assemble full generator Q and use expm
    Cap = bufCap;

    if ~isPH
        % M/M/c/N: scalar levels, (Cap+1) x (Cap+1) generator
        dim = Cap + 1;
        Q = zeros(dim, dim);
        for n = 0:Cap
            idx = n + 1;
            dep = min(n, nServers) * mu;
            arr = lambda * (n < Cap);
            if n > 0
                Q(idx, idx - 1) = dep;
            end
            if n < Cap
                Q(idx, idx + 1) = arr;
            end
            Q(idx, idx) = -(dep + arr);
        end
    else
        % M/PH/1/N: level 0 is 1-dim, levels 1..Cap are nPhases-dim
        dim = 1 + Cap * nPhases;
        Q = zeros(dim, dim);

        % Level 0
        Q(1, 1) = -lambda;
        Q(1, 2:1+nPhases) = lambda * alpha;

        for n = 1:Cap
            rows = (1 + (n-1)*nPhases + 1):(1 + n*nPhases);

            % Internal transitions (D0) and arrivals
            if n < Cap
                Q(rows, rows) = D0 - lambda * eye(nPhases);
            else
                Q(rows, rows) = D0;
            end

            % Arrivals: level n -> n+1
            if n < Cap
                next_rows = rows + nPhases;
                Q(rows, next_rows) = lambda * eye(nPhases);
            end

            % Departures: level n -> n-1
            if n == 1
                Q(rows, 1) = t_exit;
            else
                prev_rows = rows - nPhases;
                Q(rows, prev_rows) = D1;
            end
        end
    end

    % Time grid
    nTimePoints = min(101, max(11, round(T_duration * 10)));
    times = linspace(T_start, T_end, nTimePoints)';
    dt = T_duration / (nTimePoints - 1);

    % Initial distribution: empty queue
    pi0 = zeros(1, dim);
    pi0(1) = 1;

    % Compute transient distributions via iterative expm
    eQdt = expm(Q * dt);

    queue_lengths = zeros(nTimePoints, 1);
    util_values = zeros(nTimePoints, 1);
    tput_values = zeros(nTimePoints, 1);

    pi_t = pi0;
    for t_idx = 1:nTimePoints
        % Extract metrics from pi_t
        q = 0; u = 0; tput = 0;
        for n = 0:Cap
            if ~isPH
                p_n = pi_t(n + 1);
            else
                if n == 0
                    p_n = pi_t(1);
                else
                    idx_s = 1 + (n-1)*nPhases + 1;
                    idx_e = 1 + n*nPhases;
                    p_n = sum(pi_t(idx_s:idx_e));
                end
            end
            q = q + n * p_n;
            if n >= 1
                u = u + (min(n, nServers) / nServers) * p_n;
                if ~isPH
                    tput = tput + min(n, nServers) * mu * p_n;
                else
                    idx_s = 1 + (n-1)*nPhases + 1;
                    idx_e = 1 + n*nPhases;
                    tput = tput + pi_t(idx_s:idx_e) * t_exit;
                end
            end
        end
        queue_lengths(t_idx) = q;
        util_values(t_idx) = u;
        tput_values(t_idx) = tput;

        % Advance to next time point
        if t_idx < nTimePoints
            pi_t = pi_t * eQdt;
        end
    end

else
    %% Infinite capacity: use libQBD adaptive Taylor series
    proc = libqbd.QBD();

    if ~isPH
        % M/M/c: scalar QBD levels
        % Level 0: boundary
        proc.add_zero_level(-lambda, lambda);

        % Levels 1..c-1: boundary (service rate = n*mu)
        for n = 1:nServers-1
            dep = n * mu;
            proc.add_level(dep, -(lambda + dep), lambda);
        end

        % Repeating level c+: service rate = c*mu
        dep = nServers * mu;
        proc.add_final_level(dep, -(lambda + dep));
    else
        % M/PH/1: level 0 scalar, levels 1+ have nPhases phases
        proc.add_zero_level(-lambda, lambda * alpha);

        % Level 1: boundary (A_minus is m x 1, back to scalar level 0)
        A1_minus = t_exit;
        A1_0 = D0 - lambda * eye(nPhases);
        A1_plus = lambda * eye(nPhases);
        proc.add_level(A1_minus, A1_0, A1_plus);

        % Repeating level 2+: A_minus = D1 (m x m)
        proc.add_final_level(D1, D0 - lambda * eye(nPhases));
    end

    % Solve using libQBD adaptive Taylor series
    tol = options.tol;
    solver = libqbd.TaylorSeriesAdaptive(proc, {1}, tol, T_duration);

    ref_times = solver.get_reference_times() + T_start;
    ref_dists = solver.get_reference_dists();
    nTimePoints = numel(ref_times);
    times = ref_times(:);

    queue_lengths = zeros(nTimePoints, 1);
    util_values = zeros(nTimePoints, 1);
    tput_values = zeros(nTimePoints, 1);

    for t_idx = 1:nTimePoints
        dist = ref_dists{t_idx};
        nLevels = numel(dist);

        q = 0; u = 0; tput = 0;
        for n = 1:nLevels-1
            p_n = sum(dist{n + 1});
            q = q + n * p_n;
            u = u + (min(n, nServers) / nServers) * p_n;
            if ~isPH
                tput = tput + min(n, nServers) * mu * p_n;
            else
                tput = tput + dist{n + 1} * t_exit;
            end
        end
        queue_lengths(t_idx) = q;
        util_values(t_idx) = u;
        tput_values(t_idx) = tput;
    end
end

%% Package results in [metric, time] format
Qt = cell(M, K);
Ut = cell(M, K);
Tt = cell(M, K);

% Queue station
Qt{queueIdx, 1} = [queue_lengths, times];
Ut{queueIdx, 1} = [util_values, times];
Tt{queueIdx, 1} = [tput_values, times];

end
