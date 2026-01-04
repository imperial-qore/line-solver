function [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime, aoiResults] = solver_fluid_aoi(sn, options)
%SOLVER_FLUID_AOI Age of Information analysis using Markovian Fluid Queues
%
% [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime, aoiResults] = SOLVER_FLUID_AOI(sn, options)
%
% Analyzes Age of Information for single-queue open systems using the
% aoi-fluid MFQ solvers.
%
% Supported systems:
%   - Capacity 1 (bufferless): PH/PH/1/1 or PH/PH/1/1* (preemptive)
%   - Capacity 2 (single-buffer): M/PH/1/2 or M/PH/1/2* (replacement)
%
% Parameters:
%   sn (struct): Network structure from Network.getStruct()
%   options (struct): Solver options including:
%       - config.aoi_preemption: Override preemption/replacement probability (0-1)
%
% Returns:
%   QN, UN, RN, TN: Standard performance metrics
%   xvec_it, QNt, UNt, TNt, xvec_t, t: Transient outputs (for compatibility)
%   iters: Number of iterations (always 1)
%   runtime: Execution time in seconds
%   aoiResults (struct): Age of Information results:
%       .AoI_mean, .AoI_var: Mean and variance of AoI
%       .PAoI_mean, .PAoI_var: Mean and variance of Peak AoI
%       .AoI_g, .AoI_A, .AoI_h: Matrix exponential parameters for AoI CDF
%       .PAoI_g, .PAoI_A, .PAoI_h: Matrix exponential parameters for Peak AoI CDF
%       .systemType: 'bufferless' or 'singlebuffer'
%       .preemption: Preemption/replacement probability used
%
% References:
%   aoi-fluid toolbox by Ozancan Dogan, Nail Akar, Eray Unsal Atay
%   BSD 2-Clause License, 2020
%
% See also: aoi_is_aoi, aoi_extract_params, solveBufferless, solveSingleBuffer

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

runtime_start = tic;

M = sn.nstations;
K = sn.nclasses;

% Initialize output arrays
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);

% Initialize transient outputs (AoI is steady-state only)
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
xvec_it = {zeros(1, 1)};
iters = 1;

% Initialize AoI results
aoiResults = struct();
aoiResults.AoI_mean = NaN;
aoiResults.AoI_var = NaN;
aoiResults.PAoI_mean = NaN;
aoiResults.PAoI_var = NaN;
aoiResults.AoI_g = [];
aoiResults.AoI_A = [];
aoiResults.AoI_h = [];
aoiResults.PAoI_g = [];
aoiResults.PAoI_A = [];
aoiResults.PAoI_h = [];
aoiResults.systemType = '';
aoiResults.preemption = NaN;

% =========================================================================
% TOPOLOGY VALIDATION
% =========================================================================

[isAoI, aoiInfo] = aoi_is_aoi(sn);

if ~isAoI
    line_error(mfilename, 'AoI analysis requires valid AoI topology: %s', aoiInfo.errorMsg);
end

% Get station indices
sourceStation = aoiInfo.sourceStation;
queueStation = aoiInfo.queueStation;

line_debug('AoI: Source station=%d, Queue station=%d, capacity=%d', ...
    sourceStation, queueStation, aoiInfo.capacity);

% =========================================================================
% EXTRACT PARAMETERS
% =========================================================================

[aoiParams, aoiInfo] = aoi_extract_params(sn, aoiInfo, options);

aoiResults.systemType = aoiInfo.systemType;

% =========================================================================
% CALL AOI-FLUID SOLVER
% =========================================================================

if aoiInfo.capacity == 1
    % BUFFERLESS: PH/PH/1/1 or PH/PH/1/1*
    tau = aoiParams.tau;
    T = aoiParams.T;
    sigma = aoiParams.sigma;
    S = aoiParams.S;
    p = aoiParams.p;

    aoiResults.preemption = p;

    line_debug('AoI: Calling solveBufferless with p=%.2f', p);

    try
        [AoI_g, AoI_A, AoI_h, AoI_mean, AoI_var, ...
         PAoI_g, PAoI_A, PAoI_h, PAoI_mean, PAoI_var] = solveBufferless(tau, T, sigma, S, p);
    catch ME
        line_error(mfilename, 'solveBufferless failed: %s', ME.message);
    end

else
    % SINGLE-BUFFER: M/PH/1/2 or M/PH/1/2*
    lambda = aoiParams.lambda;
    sigma = aoiParams.sigma;
    S = aoiParams.S;
    r = aoiParams.r;

    aoiResults.preemption = r;

    line_debug('AoI: Calling solveSingleBuffer with r=%.2f', r);

    try
        [AoI_g, AoI_A, AoI_h, AoI_mean, AoI_var, ...
         PAoI_g, PAoI_A, PAoI_h, PAoI_mean, PAoI_var] = solveSingleBuffer(lambda, sigma, S, r);
    catch ME
        line_error(mfilename, 'solveSingleBuffer failed: %s', ME.message);
    end
end

% Store AoI results
aoiResults.AoI_mean = AoI_mean;
aoiResults.AoI_var = AoI_var;
aoiResults.PAoI_mean = PAoI_mean;
aoiResults.PAoI_var = PAoI_var;
aoiResults.AoI_g = AoI_g;
aoiResults.AoI_A = AoI_A;
aoiResults.AoI_h = AoI_h;
aoiResults.PAoI_g = PAoI_g;
aoiResults.PAoI_A = PAoI_A;
aoiResults.PAoI_h = PAoI_h;

line_debug('AoI: Mean AoI=%.4f, Var AoI=%.4f', AoI_mean, AoI_var);
line_debug('AoI: Mean PAoI=%.4f, Var PAoI=%.4f', PAoI_mean, PAoI_var);

% =========================================================================
% COMPUTE STANDARD PERFORMANCE METRICS
% =========================================================================

% Find the open class
openClasses = find(isinf(sn.njobs));
k = openClasses(1);

% Get arrival rate
lambda = sn.rates(sourceStation, k);

% Get mean service rate
if aoiInfo.capacity == 1
    % For bufferless: compute from PH representation
    % Mean service time = -sigma * S^{-1} * ones
    sigma = aoiParams.sigma;
    S = aoiParams.S;
    meanServiceTime = -sigma / S * ones(size(S, 1), 1);
    mu = 1 / meanServiceTime;
else
    % For single-buffer: compute similarly
    sigma = aoiParams.sigma;
    S = aoiParams.S;
    meanServiceTime = -sigma / S * ones(size(S, 1), 1);
    mu = 1 / meanServiceTime;
end

% Utilization
rho = lambda / mu;
UN(queueStation, k) = min(1, rho);

% Throughput (for stable system, equals arrival rate)
if rho < 1
    TN(queueStation, k) = lambda;
else
    TN(queueStation, k) = mu;  % Saturated throughput
end
TN(sourceStation, k) = TN(queueStation, k);

% Queue length (using Little's Law approximation)
% For AoI systems, we can use M/M/1 approximation for queue metrics
if rho < 1
    QN(queueStation, k) = rho / (1 - rho);
    RN(queueStation, k) = 1 / (mu - lambda);
else
    QN(queueStation, k) = Inf;
    RN(queueStation, k) = Inf;
end

% Transient outputs
QNt{queueStation, k} = [0; QN(queueStation, k)];
UNt{queueStation, k} = [0; UN(queueStation, k)];
TNt{queueStation, k} = [0; TN(queueStation, k)];

runtime = toc(runtime_start);

line_debug('AoI completed in %.4f seconds', runtime);

end
