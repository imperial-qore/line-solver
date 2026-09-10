function [QN, UN, RN, TN, pi, G] = solver_mam_bmap_map_1(bmap, service_map, options)
% SOLVER_MAM_BMAP_MAP_1 Solves a BMAP/MAP/1 queue using M/G/1 type analysis
%
% [QN, UN, RN, TN, PI, G] = SOLVER_MAM_BMAP_MAP_1(BMAP, SERVICE_MAP) solves a
% BMAP/MAP/1 queue with batch Markovian arrival process BMAP and Markovian
% service process SERVICE_MAP. Uses M/G/1 type matrix-analytic methods via
% the ETAQA algorithm.
%
% Input:
%   bmap        - BMAP object or cell array {D0, D1, D2, ..., DK} where
%                 D0: transitions without arrivals
%                 Dk: transitions triggering batch size k arrivals (k >= 1)
%   service_map - MAP object or cell array {S0, S1} where
%                 S0: transitions without service completions
%                 S1: transitions triggering service completions
%   options     - (optional) Options structure with fields:
%                 .nMoments: number of queue length moments to compute (default: 3)
%
% Output:
%   QN      - Mean queue length (1st moment)
%   UN      - Server utilization
%   RN      - Mean response time
%   TN      - Throughput (arrival rate)
%   pi      - Aggregated stationary probability [pi0, pi1, pi2+pi3+...]
%   G       - G matrix (minimal nonnegative solution)
%
% The BMAP/MAP/1 queue is modeled as an M/G/1 type Markov chain where:
% - Level corresponds to queue length
% - Phases correspond to combined (BMAP phase, MAP phase) states
%
% M/G/1 type structure (level can jump UP by any amount, DOWN by exactly 1):
%   A0 = I_ma ⊗ S1           (service completion, level -1)
%   A1 = D0 ⊗ I_ms + I_ma ⊗ S0  (phase changes, level 0)
%   A_{k+1} = D_k ⊗ I_ms     (batch arrival size k, level +k)
%
% References:
%   - Stathopoulos, Riska, Hua, Smirni: ETAQA algorithm
%   - MAMSolver: www.cs.wm.edu/MAMSolver
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

% Extract BMAP matrices
if isa(bmap, 'BMAP')
    K = bmap.getMaxBatchSize();
    D = cell(1, K+1);
    D{1} = bmap.D(0);  % D0
    for k = 1:K
        D{k+1} = bmap.D(1, k);  % Dk (batch size k)
    end
elseif iscell(bmap)
    D = bmap;
    K = length(D) - 1;  % max batch size
else
    error('BMAP must be a BMAP object or cell array {D0, D1, ..., DK}');
end

% Extract MAP service matrices
if isa(service_map, 'MAP')
    S0 = service_map.D(0);
    S1 = service_map.D(1);
elseif iscell(service_map) && length(service_map) == 2
    S0 = service_map{1};
    S1 = service_map{2};
else
    error('service_map must be a MAP object or cell array {S0, S1}');
end

% Number of phases
ma = size(D{1}, 1);  % BMAP phases
ms = size(S0, 1);    % MAP service phases
m = ma * ms;         % Combined phases per level

% Validate dimensions
for k = 1:length(D)
    if size(D{k}, 1) ~= ma || size(D{k}, 2) ~= ma
        error('All BMAP matrices must be %dx%d', ma, ma);
    end
end
if size(S0, 1) ~= ms || size(S0, 2) ~= ms || size(S1, 1) ~= ms || size(S1, 2) ~= ms
    error('MAP service matrices must be %dx%d', ms, ms);
end

%% Construct M/G/1 type matrices using Kronecker products
% A = [A0, A1, A2, ..., A_{K+1}] for levels >= 1
% A0 = I_ma ⊗ S1 (service completion, level -1)
% A1 = D0 ⊗ I_ms + I_ma ⊗ S0 (phase changes, level 0)
% A_{k+1} = D_k ⊗ I_ms for k >= 1 (batch arrivals size k, level +k)

I_ma = eye(ma);
I_ms = eye(ms);

A = zeros(m, m * (K + 2));

% A0 = I_ma ⊗ S1 (service completion)
A0 = kron(I_ma, S1);
A(:, 1:m) = A0;

% A1 = D0 ⊗ I_ms + I_ma ⊗ S0 (phase changes only)
A1 = kron(D{1}, I_ms) + kron(I_ma, S0);
A(:, m+1:2*m) = A1;

% A_{k+1} = D_k ⊗ I_ms for k >= 1 (batch arrivals)
for k = 1:K
    Ak = kron(D{k+1}, I_ms);
    A(:, (k+1)*m + 1 : (k+2)*m) = Ak;
end

% B = [B0, B1, B2, ..., BK] for level 0 (empty queue)
% At level 0, no service can occur, so we add S1 back to the diagonal block
% B0 = D0 ⊗ I_ms + I_ma ⊗ (S0 + S1)
% B_k = D_k ⊗ I_ms for k >= 1

B = zeros(m, m * (K + 1));

% B0: At level 0, no service (add S1 back as self-loop)
B0 = kron(D{1}, I_ms) + kron(I_ma, S0 + S1);
B(:, 1:m) = B0;

% B_k = D_k ⊗ I_ms for k >= 1 (batch arrivals from empty queue)
for k = 1:K
    Bk = kron(D{k+1}, I_ms);
    B(:, k*m + 1 : (k+1)*m) = Bk;
end

%% Verify stability condition
% Compute arrival rate from BMAP
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
    warning('solver_mam_bmap_map_1:Unstable', ...
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
