function [QN, UN, RN, TN, piAgg, GM] = solver_mam_qsys(arrival, service, options)
% SOLVER_MAM_QSYS Unified solver for MAP/BMAP queueing systems with single server
%
% [QN, UN, RN, TN, PI, GM] = SOLVER_MAM_QSYS(ARRIVAL, SERVICE) solves a queueing
% system with Markovian arrival and service processes. Automatically detects
% whether arrivals or service are batched and applies the appropriate
% matrix-analytic method:
%   - BMAP/MAP/1: M/G/1 type analysis (batch arrivals, single service)
%   - MAP/BMAP/1: GI/M/1 type analysis (single arrivals, batch service)
%
% Input:
%   arrival - Arrival process: MAP object, BMAP object, or cell array
%             MAP: {D0, D1} or MAP object
%             BMAP: {D0, D1, D2, ..., DK} or BMAP object
%   service - Service process: MAP object, BMAP object, or cell array
%             MAP: {S0, S1} or MAP object
%             BMAP: {S0, S1, S2, ..., SK} or BMAP object
%   options - (optional) Options structure with fields:
%             .nMoments: number of queue length moments to compute (default: 3)
%
% Output:
%   QN    - Mean queue length E[N]
%   UN    - Server utilization rho
%   RN    - Mean response time E[R]
%   TN    - Throughput (arrival rate)
%   piAgg - Aggregated stationary probabilities
%   GM    - G matrix (M/G/1 type) or R matrix (GI/M/1 type)
%
% Queue Types:
%   BMAP/MAP/1 - Batch arrivals, single service (M/G/1 type)
%     Level increases by batch size (1, 2, 3, ...), decreases by 1
%   MAP/BMAP/1 - Single arrivals, batch service (GI/M/1 type)
%     Level increases by 1, decreases by batch size (1, 2, 3, ...)
%
% References:
%   Stathopoulos, Riska, Hua, Smirni: ETAQA algorithm
%   MAMSolver: www.cs.wm.edu/MAMSolver
%
% See also: MAP, BMAP
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Input validation and setup
if nargin < 3
    options = struct();
end
if ~isfield(options, 'nMoments')
    options.nMoments = 3;
end

%% Determine arrival type (MAP or BMAP)
[arrivalMatrices, arrivalIsBatch, ma] = parseProcess(arrival, 'arrival');

%% Determine service type (MAP or BMAP)
[serviceMatrices, serviceIsBatch, ms] = parseProcess(service, 'service');

%% Route to appropriate solver based on queue type
if arrivalIsBatch && ~serviceIsBatch
    % BMAP/MAP/1 queue - M/G/1 type analysis
    [QN, UN, RN, TN, piAgg, GM] = solveBMAPMAP1(arrivalMatrices, serviceMatrices, ma, ms, options);
elseif ~arrivalIsBatch && serviceIsBatch
    % MAP/BMAP/1 queue - GI/M/1 type analysis
    [QN, UN, RN, TN, piAgg, GM] = solveMAPBMAP1(arrivalMatrices, serviceMatrices, ma, ms, options);
elseif ~arrivalIsBatch && ~serviceIsBatch
    % MAP/MAP/1 queue - treat as BMAP/MAP/1 with batch size 1
    [QN, UN, RN, TN, piAgg, GM] = solveBMAPMAP1(arrivalMatrices, serviceMatrices, ma, ms, options);
else
    % BMAP/BMAP/1 - both batched (not yet supported)
    error('BMAP/BMAP/1 queues with both batch arrivals and batch service are not yet supported');
end

end

%% Helper function to parse arrival/service process
function [matrices, isBatch, nPhases] = parseProcess(proc, procType)
    if isa(proc, 'BMAP')
        K = proc.getMaxBatchSize();
        matrices = cell(1, K + 1);
        matrices{1} = proc.D(0);  % D0
        for k = 1:K
            matrices{k+1} = proc.D(1, k);  % Dk (batch size k)
        end
        isBatch = K > 1;
        nPhases = size(matrices{1}, 1);
    elseif isa(proc, 'MAP')
        matrices = {proc.D(0), proc.D(1)};
        isBatch = false;
        nPhases = size(matrices{1}, 1);
    elseif iscell(proc)
        matrices = proc;
        isBatch = length(proc) > 2;  % More than {D0, D1} means batch process
        nPhases = size(proc{1}, 1);
    else
        error('%s must be a MAP object, BMAP object, or cell array', procType);
    end

    % Validate dimensions
    for k = 1:length(matrices)
        if size(matrices{k}, 1) ~= nPhases || size(matrices{k}, 2) ~= nPhases
            error('All %s matrices must be %dx%d', procType, nPhases, nPhases);
        end
    end
end

%% BMAP/MAP/1 solver using M/G/1 type analysis
function [QN, UN, RN, TN, pi, G] = solveBMAPMAP1(D, S, ma, ms, options)
% Solves BMAP/MAP/1 using M/G/1 type matrix-analytic methods
%
% M/G/1 type structure (level can jump UP by any amount, DOWN by exactly 1):
%   A0 = I_ma \otimes S1           (service completion, level -1)
%   A1 = D0 \otimes I_ms + I_ma \otimes S0  (phase changes, level 0)
%   A_{k+1} = D_k \otimes I_ms     (batch arrival size k, level +k)

K = length(D) - 1;  % max batch size for arrivals
m = ma * ms;        % Combined phases per level

% Extract MAP service matrices
S0 = S{1};
S1 = S{2};

%% Construct M/G/1 type matrices using Kronecker products
I_ma = eye(ma);
I_ms = eye(ms);

A = zeros(m, m * (K + 2));

% A0 = I_ma \otimes S1 (service completion)
A0 = kron(I_ma, S1);
A(:, 1:m) = A0;

% A1 = D0 \otimes I_ms + I_ma \otimes S0 (phase changes only)
A1 = kron(D{1}, I_ms) + kron(I_ma, S0);
A(:, m+1:2*m) = A1;

% A_{k+1} = D_k \otimes I_ms for k >= 1 (batch arrivals)
for k = 1:K
    Ak = kron(D{k+1}, I_ms);
    A(:, (k+1)*m + 1 : (k+2)*m) = Ak;
end

% B = [B0, B1, B2, ..., BK] for level 0 (empty queue)
B = zeros(m, m * (K + 1));

% B0: At level 0, no service (add S1 back as self-loop)
B0 = kron(D{1}, I_ms) + kron(I_ma, S0 + S1);
B(:, 1:m) = B0;

% B_k = D_k \otimes I_ms for k >= 1 (batch arrivals from empty queue)
for k = 1:K
    Bk = kron(D{k+1}, I_ms);
    B(:, k*m + 1 : (k+1)*m) = Bk;
end

%% Verify stability condition
D0 = D{1};
D1_total = zeros(ma, ma);
for k = 1:K
    D1_total = D1_total + D{k+1};
end

% Stationary distribution of BMAP
pi_bmap = map_prob({D0, D1_total});
e_ma = ones(ma, 1);

% Total customer arrival rate: sum_k (k * π * Dk * e)
lambda_total = 0;
for k = 1:K
    rate_k = pi_bmap * D{k+1} * e_ma;
    lambda_total = lambda_total + k * rate_k;
end

% Compute service rate from MAP
pi_map = map_prob({S0, S1});
e_ms = ones(ms, 1);
mu = pi_map * S1 * e_ms;  % Service rate

% Utilization
rho = lambda_total / mu;

if rho >= 1
    warning('solver_mam_qsys:Unstable', ...
        'System is unstable (rho = %.4f >= 1). Results may be invalid.', rho);
end

%% Compute G matrix using MG1_G_ETAQA
G = MG1_G_ETAQA(A);

%% Compute stationary probabilities using MG1_pi_ETAQA
pi = MG1_pi_ETAQA(B, A, G);

%% Compute queue length moments using MG1_qlen_ETAQA
qlen_moments = zeros(1, options.nMoments);
for n = 1:options.nMoments
    qlen_moments(n) = MG1_qlen_ETAQA(B, A, pi, n);
end

%% Compute performance metrics
QN = qlen_moments(1);  % Mean queue length (E[N])
UN = rho;              % Utilization
TN = lambda_total;     % Throughput (total arrival rate)
RN = QN / TN;          % Mean response time (Little's law: E[R] = E[N] / λ)

end

%% MAP/BMAP/1 solver using GI/M/1 type analysis
function [QN, UN, RN, TN, piAgg, R] = solveMAPBMAP1(C, D, ma, ms, options)
% Solves MAP/BMAP/1 using GI/M/1 type matrix-analytic methods
%
% GI/M/1 type structure (level jumps UP by 1, DOWN by any amount):
%   A0 = C1 \otimes I_ms           (MAP arrival, level +1)
%   A1 = C0 \otimes I_ms + I_ma \otimes D0  (phase changes, level 0)
%   A_{k+1} = I_ma \otimes D_k     (batch size k service, level -k)

C0 = C{1};
C1 = C{2};
K = length(D) - 1;  % max batch size for service
m = ma * ms;        % Combined phases per level

%% Validate inputs
% Check MAP is valid generator
rowSumsC = sum(C0 + C1, 2);
if max(abs(rowSumsC)) > 1e-10
    error('MAP matrices C0 + C1 must have zero row sums');
end

% Check BMAP is valid generator
DTotal = D{1};
for k = 1:K
    DTotal = DTotal + D{k+1};
end
rowSumsD = sum(DTotal, 2);
if max(abs(rowSumsD)) > 1e-10
    error('BMAP matrices D0 + D1 + ... + DK must have zero row sums');
end

%% Compute arrival and service rates
% Arrival rate from MAP
piC = map_prob({C0, C1});
lambdaArr = piC * C1 * ones(ma, 1);

% Service rate from BMAP (total customers served per unit time)
piD = map_prob({D{1}, DTotal - D{1}});  % D0, sum of D1...DK
muTotal = 0;
for k = 1:K
    muTotal = muTotal + k * (piD * D{k+1} * ones(ms, 1));
end

% Utilization
rho = lambdaArr / muTotal;

if rho >= 1.0
    warning('solver_mam_qsys:Unstable', ...
        'System is unstable (rho = %.4f >= 1). Results may be invalid.', rho);
end

%% Construct A matrix for GI/M/1-type (repeating part)
A = zeros(m * (K + 2), m);

% A0: MAP arrival (level +1)
A0 = kron(C1, eye(ms));
A(1:m, :) = A0;

% A1: Phase changes only (level 0)
A1 = kron(C0, eye(ms)) + kron(eye(ma), D{1});
A(m+1:2*m, :) = A1;

% A_{k+1}: Batch service of size k (level -k)
for k = 1:K
    Ak = kron(eye(ma), D{k+1});
    A((k+1)*m+1:(k+2)*m, :) = Ak;
end

%% Construct B matrix for boundary levels
B = zeros(m * (K + 1) + m, m);

% B1 (level 0 -> level 0): Phase changes only, no service from empty
B1 = kron(C0, eye(ms)) + kron(eye(ma), D{1});
% Add service transitions that would go to negative levels back to level 0
for k = 1:K
    B1 = B1 + kron(eye(ma), D{k+1});
end
B(1:m, :) = B1;

% B2, B3, ..., B_{K+1}: Transitions to level 0 from levels 1, 2, ..., K
for j = 1:K
    Bj = zeros(m, m);
    for k = j:K
        Bj = Bj + kron(eye(ma), D{k+1});
    end
    B(m + (j-1)*m + 1 : m + j*m, :) = Bj;
end

%% Compute R matrix using GI/M/1 ETAQA
R = GIM1_R_ETAQA(A);

%% Compute stationary probabilities using GI/M/1 ETAQA
piAgg = GIM1_pi_ETAQA(B, A, R, 'Boundary', A0);

%% Compute performance metrics
% Mean queue length using ETAQA moments
QN = GIM1_qlen_ETAQA(B, A, R, piAgg, 1, 'Boundary', A0);

% Throughput equals arrival rate (in steady state)
TN = lambdaArr;

% Utilization
UN = rho;

% Mean response time via Little's Law
RN = QN / TN;

end
