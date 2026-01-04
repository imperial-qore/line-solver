function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_diffusion(sn, options)
% SOLVER_FLUID_DIFFUSION Diffusion approximation for closed multiclass BCMP networks
%
% [QN,UN,RN,TN,XVEC_IT,QNT,UNT,TNT,XVEC_T,T,ITERS,RUNTIME] = SOLVER_FLUID_DIFFUSION(SN, OPTIONS)
%
% Implements diffusion approximation using the Euler-Maruyama method for SDEs.
% This method models queue lengths as continuous variables with Brownian noise,
% providing stochastic extensions to deterministic fluid approximations.
%
% Supported models:
%   - Closed multiclass queueing networks only
%   - Scheduling: PS, FCFS, INF (delay stations)
%   - Single-server or infinite-server stations
%
% Options:
%   - iter_max: Number of simulation steps (default: 10000)
%   - timestep: Time step for Euler-Maruyama integration (default: 0.01)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

runtime_start = tic;

M = sn.nstations;    % number of stations
K = sn.nclasses;     % number of classes
N_vec = sn.njobs';   % population per class (column vector)

%% Model validation - diffusion only supports closed networks
if any(isinf(N_vec))
    line_error(mfilename, 'Diffusion method only supports closed queueing networks (no open classes).');
end

% Check for sources or sinks (not allowed in closed networks for diffusion)
for ist = 1:M
    if sn.sched(ist) == SchedStrategy.EXT
        line_error(mfilename, 'Diffusion method does not support Source nodes. Use closed networks only.');
    end
end

% Check scheduling disciplines - only PS, FCFS, INF supported
for ist = 1:M
    switch sn.sched(ist)
        case {SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.INF, SchedStrategy.SIRO}
            % supported
        otherwise
            line_error(mfilename, sprintf('Diffusion method does not support scheduling strategy %s at station %d.', ...
                SchedStrategy.toText(sn.sched(ist)), ist));
    end
end

% Check for single-server or infinite-server only
for ist = 1:M
    if sn.nservers(ist) > 1 && ~isinf(sn.nservers(ist))
        line_error(mfilename, sprintf('Diffusion method only supports single-server (c=1) or infinite-server stations. Station %d has %d servers.', ...
            ist, sn.nservers(ist)));
    end
end

%% Extract diffusion-specific options from existing solver options
% Use iter_max for number of simulation steps
if isfield(options, 'iter_max') && ~isempty(options.iter_max)
    steps = options.iter_max;
else
    steps = 10000;  % default number of steps
end

% Use timestep for dt
if isfield(options, 'timestep') && ~isempty(options.timestep)
    dt = options.timestep;
else
    dt = 0.01;  % default time step
end

% Total simulation time is derived from steps and dt
T = steps * dt;

%% Build inverse service rate matrix mu_inv (M x K)
% mu_inv(k,r) = average service time for class r at station k
mu_inv = zeros(M, K);
for ist = 1:M
    for r = 1:K
        if sn.rates(ist, r) > 0 && ~isnan(sn.rates(ist, r))
            mu_inv(ist, r) = 1 / sn.rates(ist, r);
        else
            % If no service for this class at this station, set to Inf
            % (effectively no departures)
            mu_inv(ist, r) = Inf;
        end
    end
end

%% Build routing matrix P ((M*K) x (M*K))
% P(i,j) is probability of routing from flattened state i to flattened state j
% Flattening: index = (station-1)*K + class for station-major ordering
% However, bcmp_diffusion uses sub2ind([K R], k, r) which is column-major [queue, class]
% So we use sub2ind([M K], ist, r) for consistency

% Extract station-to-station routing from the full routing table
P_full = sn.rt;  % Full routing matrix (stateful nodes, class-blocked)

% Build station indices mapping
station_indices = [];
for ist = 1:M
    isf = sn.stationToStateful(ist);
    for r = 1:K
        station_indices = [station_indices, (isf-1)*K + r];
    end
end

% Extract the station-level routing matrix
P_stations = P_full(station_indices, station_indices);

% Reshape to (M*K) x (M*K) in the format expected by diffusion
% The routing matrix needs to be in [station, class] flattened order
P = zeros(M*K, M*K);
for ist_from = 1:M
    for r_from = 1:K
        idx_from = sub2ind([M K], ist_from, r_from);
        for ist_to = 1:M
            for r_to = 1:K
                idx_to = sub2ind([M K], ist_to, r_to);
                % Map from station-blocked format to sub2ind format
                src_idx = (ist_from-1)*K + r_from;
                dst_idx = (ist_to-1)*K + r_to;
                P(idx_from, idx_to) = P_stations(src_idx, dst_idx);
            end
        end
    end
end

%% Run diffusion approximation
line_debug('Diffusion solver starting: steps=%d, dt=%g, nstations=%d, nclasses=%d', steps, dt, M, K);

[Xavg, Xt] = solver_fluid_diffusion_core(mu_inv, P, N_vec, T, dt);

runtime = toc(runtime_start);
iters = 1;  % Diffusion runs a single trajectory

%% Extract results
steps = size(Xt, 3);
t = (0:steps-1)' * dt;

% Queue lengths (steady-state average)
QN = Xavg;

% Transient queue lengths
QNt = cell(M, K);
for ist = 1:M
    for r = 1:K
        QNt{ist, r} = squeeze(Xt(ist, r, :));
    end
end

% Throughputs: TN = QN * service_rate (for PS/INF) using Little's law approximation
TN = zeros(M, K);
for ist = 1:M
    for r = 1:K
        if mu_inv(ist, r) > 0 && ~isinf(mu_inv(ist, r))
            % For infinite servers, throughput = queue_length * service_rate
            % For single server PS, this is an approximation
            if isinf(sn.nservers(ist))
                TN(ist, r) = QN(ist, r) / mu_inv(ist, r);
            else
                % For single server, throughput limited by server capacity
                TN(ist, r) = min(QN(ist, r), 1) / mu_inv(ist, r);
            end
        end
    end
end

% Transient throughputs
TNt = cell(M, K);
for ist = 1:M
    for r = 1:K
        if mu_inv(ist, r) > 0 && ~isinf(mu_inv(ist, r))
            if isinf(sn.nservers(ist))
                TNt{ist, r} = QNt{ist, r} / mu_inv(ist, r);
            else
                TNt{ist, r} = min(QNt{ist, r}, 1) / mu_inv(ist, r);
            end
        else
            TNt{ist, r} = zeros(steps, 1);
        end
    end
end

% Utilizations
UN = zeros(M, K);
UNt = cell(M, K);
for ist = 1:M
    S = sn.nservers(ist);
    for r = 1:K
        if isinf(S)
            % Infinite server: utilization = queue length
            UN(ist, r) = QN(ist, r);
            UNt{ist, r} = QNt{ist, r};
        else
            % Single server: utilization = min(queue_length, 1)
            UN(ist, r) = min(QN(ist, r) / S, 1);
            UNt{ist, r} = min(QNt{ist, r} / S, 1);
        end
    end
end

% Response times using Little's law: RN = QN / TN
RN = zeros(M, K);
for ist = 1:M
    for r = 1:K
        if TN(ist, r) > 0
            RN(ist, r) = QN(ist, r) / TN(ist, r);
        end
    end
end

% State vectors for compatibility
xvec_it = {Xavg(:)};
xvec_t = Xt;

line_debug('Diffusion solver completed in %g seconds', runtime);

end

%% Core diffusion algorithm
function [Xavg, Xt] = solver_fluid_diffusion_core(mu_inv, P, N_vec, T, dt)
% SOLVER_FLUID_DIFFUSION_CORE Core diffusion approximation for BCMP/Kelly networks
%
% Inputs:
%   mu_inv : (M x K) matrix of inverse service rates (average service time)
%   P      : ((M*K) x (M*K)) routing matrix (flattened, column-major)
%   N_vec  : (K x 1) vector of total jobs per class
%   T      : total simulation time
%   dt     : time step
%
% Output:
%   Xavg   : (M x K) matrix of time-averaged queue lengths
%   Xt     : (M x K x steps) tensor of queue lengths over time

[M, K] = size(mu_inv);
steps = floor(T/dt);

% Initialize storage
Xt = zeros(M, K, steps);

% Initial state: distribute jobs evenly across stations for each class
for r = 1:K
    Xt(:, r, 1) = N_vec(r) / M;
end

% Initialize running sum for time-average
Xavg = Xt(:, :, 1) / steps;

% Drift function handle
drift_fn = @(x) diffusion_drift(x, mu_inv, P, M, K);

% Euler-Maruyama simulation
for step = 2:steps
    x_prev = Xt(:, :, step-1);

    % Brownian motion increment
    dW = sqrt(dt) * randn(M, K);

    % Drift term
    drift = drift_fn(x_prev);

    % Euler-Maruyama update
    x_new = x_prev + drift * dt + dW;

    % Reflect at zero (queue lengths cannot be negative)
    x_new = max(x_new, 0);

    % Project back to fixed population per class (closed network constraint)
    for r = 1:K
        total_r = sum(x_new(:, r));
        if total_r > 0
            x_new(:, r) = x_new(:, r) * N_vec(r) / total_r;
        else
            % Re-distribute evenly if all became zero
            x_new(:, r) = N_vec(r) / M;
        end
    end

    % Store state
    Xt(:, :, step) = x_new;

    % Accumulate for time-average
    Xavg = Xavg + x_new / steps;
end

end

%% Drift function
function d = diffusion_drift(x, mu_inv, P, M, K)
% DIFFUSION_DRIFT Calculate net flow (drift) for each (station, class) pair
%
% Inputs:
%   x       : M x K matrix of current queue lengths
%   mu_inv  : M x K matrix of average service times
%   P       : (M*K x M*K) routing matrix
%   M       : Number of stations
%   K       : Number of classes
%
% Output:
%   d       : M x K matrix of drift values

d = zeros(M, K);

for ist = 1:M
    for r = 1:K
        idx_out = sub2ind([M K], ist, r);

        % Outflow: departure rate from (station, class)
        if mu_inv(ist, r) > 0 && ~isinf(mu_inv(ist, r))
            out_flow = x(ist, r) / mu_inv(ist, r);
        else
            out_flow = 0;
        end

        % Inflow: sum of arrivals from all other (station, class) pairs
        inflow = 0;
        for j = 1:(M*K)
            [ist_j, r_j] = ind2sub([M K], j);
            if mu_inv(ist_j, r_j) > 0 && ~isinf(mu_inv(ist_j, r_j))
                inflow = inflow + (x(ist_j, r_j) / mu_inv(ist_j, r_j)) * P(j, idx_out);
            end
        end

        % Net drift
        d(ist, r) = inflow - out_flow;
    end
end

end
