function result = qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, varargin)
% QSYS_BMAPPHNN_RETRIAL Analyzes a BMAP/PH/N/N bufferless retrial queue.
%
% RESULT = QSYS_BMAPPHNN_RETRIAL(D, BETA, S, N, ALPHA, GAMMA, P, R) analyzes
% a BMAP/PH/N/N bufferless retrial queueing system with admission control.
%
% This implements the algorithm from:
% Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System with Flexible
% Retrials Admission Control", Mathematics 2025, 13(9), 1434.
%
% Inputs:
%   D     - Cell array {D0, D1, ..., DK} of BMAP matrices
%           D0: hidden transition matrix (V x V)
%           D1, ..., DK: arrival matrices for batches of size 1, ..., K
%   BETA  - PH service initial probability vector (1 x M)
%   S     - PH service subgenerator matrix (M x M)
%   N     - Number of servers (also capacity, hence bufferless)
%   ALPHA - Retrial rate per customer in orbit
%   GAMMA - Impatience (abandonment) rate per customer in orbit
%   P     - Probability of batch rejection when not enough servers
%   R     - Admission threshold (scalar or 1 x V vector per BMAP state)
%           When n > R(nu), arriving customers go to orbit
%
% Optional parameters:
%   'MaxLevel'  - Maximum orbit level for truncation (default: auto)
%   'Tolerance' - Convergence tolerance (default: 1e-10)
%   'Verbose'   - Print progress messages (default: false)
%
% Returns a struct with fields:
%   L_orbit        - Mean number of customers in orbit
%   N_server       - Mean number of busy servers
%   L_system       - Mean number in system (orbit + servers)
%   Utilization    - Server utilization (N_server / N)
%   Throughput     - System throughput
%   P_idle         - Probability all servers are idle
%   P_empty_orbit  - Probability orbit is empty
%   P_empty_system - Probability system is empty (idle and empty orbit)
%   pi             - Stationary distribution (levels x Vd)
%   truncLevel     - Truncation level used
%   analyzer       - Name of analyzer used
%
% Example:
%   % M/M/3/3 retrial queue (exponential arrivals and service)
%   D = {-2.0, 2.0};  % Exp(2) arrivals
%   beta = 1;
%   S = -1;           % Exp(1) service
%   N = 3;
%   alpha = 0.5;      % Retrial rate
%   gamma = 0;        % No impatience
%   p = 0;            % No batch rejection
%   R = 2;            % Admission threshold
%   result = qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R);
%
% See also qsys_mapph1, qsys_is_retrial
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Parse optional arguments
parser = inputParser;
addParameter(parser, 'MaxLevel', [], @isnumeric);
addParameter(parser, 'Tolerance', 1e-10, @isnumeric);
addParameter(parser, 'Verbose', false, @islogical);
parse(parser, varargin{:});

maxLevelParam = parser.Results.MaxLevel;
tol = parser.Results.Tolerance;
verbose = parser.Results.Verbose;

%% Validate and process inputs

% BMAP parameters
K = length(D) - 1;  % Maximum batch size
V = size(D{1}, 1);  % Number of BMAP states

% Compute generator of fundamental process: D^(1) = sum(D_k)
D1_gen = zeros(size(D{1}));
for k = 1:length(D)
    D1_gen = D1_gen + D{k};
end

% Stationary distribution of fundamental process
theta = computeStationaryVector(D1_gen);

% Mean arrival rate: lambda = theta * sum(k * D_k) * e
sumKDk = zeros(size(D{1}));
for k = 2:length(D)
    sumKDk = sumKDk + (k-1) * D{k};
end
lambda = theta * sumKDk * ones(V, 1);

% PH service parameters
beta = beta(:)';
M = size(S, 1);
S0 = -S * ones(M, 1);
b1 = beta * (-S \ ones(M, 1));  % Mean service time

% Handle R parameter (threshold)
if isscalar(R)
    R = R * ones(1, V);
else
    R = R(:)';
end

% Compute T_n values: T_n = C(n+M-1, M-1) = number of service states with n busy
T = zeros(1, N+1);
for n = 0:N
    T(n+1) = nchoosek(n + M - 1, M - 1);
end
d = sum(T);  % Total dimension per BMAP state

% Build state mapping
stateMap = buildStateMap(N, M, T);

%% Determine truncation level
rho = lambda * b1 / N;  % Offered load

if verbose
    fprintf('Solving BMAP/PH/N/N retrial queue...\n');
    fprintf('  V=%d, M=%d, N=%d, K=%d\n', V, M, N, K);
    fprintf('  d=%d, block size Vd=%d\n', d, V*d);
    fprintf('  lambda=%.4f, mu=%.4f\n', lambda, 1/b1);
    fprintf('  Offered load rho=%.4f\n', rho);
end

if isempty(maxLevelParam)
    truncLevel = max(100, ceil(50 / (1 - min(rho, 0.99))));
else
    truncLevel = maxLevelParam;
end

if verbose
    fprintf('  Truncation level: %d\n', truncLevel);
end

%% Build and solve the system
Vd = V * d;
totalDim = (truncLevel + 1) * Vd;

if verbose
    fprintf('Total matrix dimension: %d x %d\n', totalDim, totalDim);
end

% Use sparse for large systems
if totalDim > 5000
    if verbose
        fprintf('Using sparse representation\n');
    end
    Q = sparse(totalDim, totalDim);
else
    Q = zeros(totalDim, totalDim);
end

% Context structure for helper functions
ctx = struct('D', {D}, 'beta', beta, 'S', S, 'S0', S0, 'M', M, 'N', N, ...
    'V', V, 'K', K, 'd', d, 'T', T, 'R', R, 'alpha', alpha, ...
    'gamma', gamma, 'p', p, 'stateMap', {stateMap});

% Fill generator blocks
for i = 0:truncLevel
    if verbose && mod(i, 20) == 0
        fprintf('  Level %d/%d\n', i, truncLevel);
    end

    rowStart = i*Vd + 1;
    rowEnd = (i+1)*Vd;

    for j = max(0, i-1):min(truncLevel, i+K)
        colStart = j*Vd + 1;
        colEnd = (j+1)*Vd;

        Qij = buildGeneratorLevel(ctx, i, j);
        Q(rowStart:rowEnd, colStart:colEnd) = Qij;
    end
end

% Ensure rows sum to zero
for i = 1:totalDim
    Q(i,i) = Q(i,i) - sum(Q(i,:));
end

% Solve pi * Q = 0, pi * e = 1
if verbose
    fprintf('Solving linear system...\n');
end

Q(:, end) = ones(totalDim, 1);
b = zeros(1, totalDim);
b(end) = 1;

if issparse(Q)
    pi = (Q' \ b')';
else
    pi = b / Q;
end

% Reshape to level structure
pi = reshape(pi, Vd, truncLevel + 1)';

% Handle numerical issues
if any(pi(:) < -1e-8)
    line_warning(mfilename, 'Negative probabilities detected, clipping to zero');
    pi(pi < 0) = 0;
end
pi = pi / sum(pi(:));  % Renormalize

%% Compute performance measures
maxLevel = size(pi, 1) - 1;

% Mean number in orbit
L_orbit = 0;
for i = 1:maxLevel
    L_orbit = L_orbit + i * sum(pi(i+1, :));
end

% Mean number of busy servers
N_server = 0;
for i = 0:maxLevel
    piLevel = pi(i+1, :);
    for nu = 1:V
        for n = 0:N
            offset = (nu-1)*d + getBlockOffset(T, n);
            for t = 1:T(n+1)
                idx = offset + t - 1;
                if idx <= Vd
                    N_server = N_server + n * piLevel(idx);
                end
            end
        end
    end
end

% Probability all servers idle
P_idle = 0;
for i = 0:maxLevel
    piLevel = pi(i+1, :);
    for nu = 1:V
        offset = (nu-1)*d + 1;  % n=0
        P_idle = P_idle + piLevel(offset);
    end
end

% Probability orbit empty
P_empty_orbit = sum(pi(1, :));

% Probability system empty
P_empty = 0;
piLevel = pi(1, :);
for nu = 1:V
    offset = (nu-1)*d + 1;
    P_empty = P_empty + piLevel(offset);
end

%% Build result struct
result = struct();
result.L_orbit = L_orbit;
result.N_server = N_server;
result.L_system = L_orbit + N_server;
result.Utilization = N_server / N;
result.Throughput = N_server / b1;
result.P_idle = P_idle;
result.P_empty_orbit = P_empty_orbit;
result.P_empty_system = P_empty;
result.pi = pi;
result.truncLevel = truncLevel;
result.analyzer = 'LINE:qsys_bmapphnn_retrial';

if verbose
    fprintf('Solution complete.\n');
end

end

%% ========== Helper Functions ==========

function theta = computeStationaryVector(Q)
% Solve theta * Q = 0, theta * e = 1
n = size(Q, 1);
A = Q';
A(end, :) = ones(1, n);
b = zeros(n, 1);
b(end) = 1;
theta = (A \ b)';
end

function stateMap = buildStateMap(N, M, T)
% Build mapping from (n, service_state_vector) to linear index
stateMap = cell(N + 1, 1);
for n = 0:N
    stateMap{n+1} = generateCompositions(n, M);
end
end

function comps = generateCompositions(n, M)
% Generate all weak compositions of n into M parts (reverse lexicographic)
if M == 1
    comps = n;
    return;
end
numComps = nchoosek(n + M - 1, M - 1);
comps = zeros(numComps, M);
idx = 1;
for m1 = n:-1:0
    subComps = generateCompositions(n - m1, M - 1);
    numSub = size(subComps, 1);
    comps(idx:idx+numSub-1, 1) = m1;
    comps(idx:idx+numSub-1, 2:end) = subComps;
    idx = idx + numSub;
end
end

function offset = getBlockOffset(T, n)
% Get starting index (1-based) for states with n busy servers
if n == 0
    offset = 1;
else
    offset = sum(T(1:n)) + 1;
end
end

function L = computeL(ctx, n)
% Matrix L_n: service completion transitions (T_n x T_{n-1})
if n == 0
    L = [];
    return;
end
L = zeros(ctx.T(n+1), ctx.T(n));
compsN = ctx.stateMap{n+1};
compsNm1 = ctx.stateMap{n};

for i = 1:size(compsN, 1)
    m = compsN(i, :);
    for l = 1:ctx.M
        if m(l) > 0
            mPrime = m;
            mPrime(l) = mPrime(l) - 1;
            for j = 1:size(compsNm1, 1)
                if all(compsNm1(j,:) == mPrime)
                    L(i, j) = L(i, j) + m(l) * ctx.S0(l);
                    break;
                end
            end
        end
    end
end
end

function A = computeA(ctx, n)
% Matrix A_n: phase change transitions (T_n x T_n)
if n == 0
    A = 0;
    return;
end
A = zeros(ctx.T(n+1), ctx.T(n+1));
comps = ctx.stateMap{n+1};

for i = 1:size(comps, 1)
    m = comps(i, :);
    for l = 1:ctx.M
        if m(l) > 0
            for lPrime = 1:ctx.M
                if lPrime ~= l && ctx.S(l, lPrime) > 0
                    mPrime = m;
                    mPrime(l) = mPrime(l) - 1;
                    mPrime(lPrime) = mPrime(lPrime) + 1;
                    for j = 1:size(comps, 1)
                        if all(comps(j,:) == mPrime)
                            A(i, j) = A(i, j) + m(l) * ctx.S(l, lPrime);
                            break;
                        end
                    end
                end
            end
        end
    end
end
end

function P = computeP(ctx, n)
% Matrix P_n: new arrival transitions (T_n x T_{n+1})
if n >= ctx.N
    P = [];
    return;
end
P = zeros(ctx.T(n+1), ctx.T(n+2));
compsN = ctx.stateMap{n+1};
compsNp1 = ctx.stateMap{n+2};

for i = 1:size(compsN, 1)
    m = compsN(i, :);
    for l = 1:ctx.M
        if ctx.beta(l) > 0
            mPrime = m;
            mPrime(l) = mPrime(l) + 1;
            for j = 1:size(compsNp1, 1)
                if all(compsNp1(j,:) == mPrime)
                    P(i, j) = P(i, j) + ctx.beta(l);
                    break;
                end
            end
        end
    end
end
end

function Delta = computeDelta(ctx, n)
% Diagonal matrix Delta_n: exit rates (T_n x T_n)
if n == 0
    Delta = 0;
    return;
end
comps = ctx.stateMap{n+1};
diagVals = zeros(ctx.T(n+1), 1);
for i = 1:size(comps, 1)
    m = comps(i, :);
    total = 0;
    for l = 1:ctx.M
        total = total + m(l) * (-ctx.S(l, l));
    end
    diagVals(i) = total;
end
Delta = diag(diagVals);
end

function Gamma = computeGamma(ctx, nu)
% Diagonal matrix Gamma^(nu): 0 for n <= R_nu, 1 for n > R_nu
diagVals = zeros(ctx.d, 1);
offset = 1;
for n = 0:ctx.N
    if n > ctx.R(nu)
        diagVals(offset:offset + ctx.T(n+1) - 1) = 1;
    end
    offset = offset + ctx.T(n+1);
end
Gamma = diag(diagVals);
end

function G = computeG_nn(ctx, n, nu, nuPrime)
% G_{n,n}^{(nu,nu')} matrix for batch losses
if n <= ctx.N - ctx.K
    G = zeros(ctx.T(n+1));
else
    total = 0;
    for k = (ctx.N - n + 1):ctx.K
        if k >= 1 && k <= ctx.K
            total = total + ctx.D{k+1}(nu, nuPrime);
        end
    end
    G = ctx.p * total * eye(ctx.T(n+1));
end
end

function B = computeB(ctx, nu)
% Block matrix B^(nu) of size d x d
B = zeros(ctx.d, ctx.d);

% Precompute matrices
L = cell(ctx.N + 1, 1);
A = cell(ctx.N + 1, 1);
P = cell(ctx.N + 1, 1);
Delta = cell(ctx.N + 1, 1);

for n = 0:ctx.N
    L{n+1} = computeL(ctx, n);
    A{n+1} = computeA(ctx, n);
    P{n+1} = computeP(ctx, n);
    Delta{n+1} = computeDelta(ctx, n);
end

for n = 0:ctx.N
    rowStart = getBlockOffset(ctx.T, n);
    rowEnd = rowStart + ctx.T(n+1) - 1;

    % Diagonal block
    G_nn = computeG_nn(ctx, n, nu, nu);
    if n == 0
        B(rowStart, rowStart) = G_nn;
    else
        B(rowStart:rowEnd, rowStart:rowEnd) = A{n+1} + Delta{n+1} + G_nn;
    end

    % Subdiagonal block
    if n >= 1
        colStart = getBlockOffset(ctx.T, n-1);
        colEnd = colStart + ctx.T(n) - 1;
        B(rowStart:rowEnd, colStart:colEnd) = L{n+1};
    end

    % Superdiagonal blocks
    for k = 1:ctx.K
        if n + k <= ctx.N
            colStart = getBlockOffset(ctx.T, n+k);
            colEnd = colStart + ctx.T(n+k+1) - 1;
            D_k_nu_nu = ctx.D{k+1}(nu, nu);

            Pprod = eye(ctx.T(n+1));
            for j = n:n+k-1
                if j < ctx.N
                    Pprod = Pprod * P{j+1};
                end
            end
            B(rowStart:rowEnd, colStart:colEnd) = D_k_nu_nu * Pprod;
        end
    end
end
end

function Bbar = computeBbar(ctx, nu)
% B_bar^(nu) matrix for successful retrials
Bbar = zeros(ctx.d, ctx.d);
for n = 0:min(ctx.R(nu), ctx.N - 1)
    rowStart = getBlockOffset(ctx.T, n);
    rowEnd = rowStart + ctx.T(n+1) - 1;
    colStart = getBlockOffset(ctx.T, n+1);
    colEnd = colStart + ctx.T(n+2) - 1;
    P_n = computeP(ctx, n);
    Bbar(rowStart:rowEnd, colStart:colEnd) = P_n;
end
end

function Btilde = computeBtilde(ctx, nu, nuPrime)
% B_tilde^(nu, nu') for BMAP state transitions
Btilde = zeros(ctx.d, ctx.d);

P = cell(ctx.N + 1, 1);
for n = 0:ctx.N
    P{n+1} = computeP(ctx, n);
end

for n = 0:ctx.N
    rowStart = getBlockOffset(ctx.T, n);
    rowEnd = rowStart + ctx.T(n+1) - 1;

    % Diagonal block
    G_nn = computeG_nn(ctx, n, nu, nuPrime);
    Btilde(rowStart:rowEnd, rowStart:rowEnd) = G_nn;

    % Superdiagonal blocks
    for k = 1:ctx.K
        if n + k <= ctx.N
            colStart = getBlockOffset(ctx.T, n+k);
            colEnd = colStart + ctx.T(n+k+1) - 1;
            D_k_nu_nuPrime = ctx.D{k+1}(nu, nuPrime);

            Pprod = eye(ctx.T(n+1));
            for j = n:n+k-1
                if j < ctx.N
                    Pprod = Pprod * P{j+1};
                end
            end
            Btilde(rowStart:rowEnd, colStart:colEnd) = D_k_nu_nuPrime * Pprod;
        end
    end
end
end

function C = computeC(ctx, n, k, nu, nuPrime)
% C_{n,k}^(nu, nu') for partial batch admission to orbit
if n < ctx.N - ctx.K + k
    C = zeros(ctx.T(n+1), 1);
    if ctx.T(n+1) > 1
        C = zeros(ctx.T(n+1), 1);
    end
elseif n < ctx.N
    batchSize = ctx.N - n + k;
    if batchSize >= 1 && batchSize <= ctx.K
        D_batch = ctx.D{batchSize + 1}(nu, nuPrime);

        P = cell(ctx.N + 1, 1);
        for nn = 0:ctx.N
            P{nn+1} = computeP(ctx, nn);
        end

        Pprod = eye(ctx.T(n+1));
        for j = n:ctx.N-1
            Pprod = Pprod * P{j+1};
        end
        C = (1 - ctx.p) * D_batch * Pprod;
    else
        C = zeros(ctx.T(n+1), ctx.T(ctx.N+1));
    end
else  % n == N
    if k >= 1 && k <= ctx.K
        D_k = ctx.D{k+1}(nu, nuPrime);
        C = (1 - ctx.p) * D_k * eye(ctx.T(ctx.N+1));
    else
        C = zeros(ctx.T(ctx.N+1));
    end
end
end

function Q = buildGeneratorLevel(ctx, i, j)
% Build generator block Q_{i,j}
Vd = ctx.V * ctx.d;
Q = zeros(Vd, Vd);

if j < max(0, i-1) || j > i + ctx.K
    return;
end

% Precompute matrices
B = cell(ctx.V, 1);
Bbar = cell(ctx.V, 1);
Gamma = cell(ctx.V, 1);

for nu = 1:ctx.V
    B{nu} = computeB(ctx, nu);
    Bbar{nu} = computeBbar(ctx, nu);
    Gamma{nu} = computeGamma(ctx, nu);
end

if i == j  % Diagonal block
    for nu = 1:ctx.V
        rowStart = (nu-1)*ctx.d + 1;
        rowEnd = nu*ctx.d;

        for nuPrime = 1:ctx.V
            colStart = (nuPrime-1)*ctx.d + 1;
            colEnd = nuPrime*ctx.d;

            if nu == nuPrime
                D0_nu_nu = ctx.D{1}(nu, nu);
                block = D0_nu_nu * eye(ctx.d) + B{nu} ...
                    - i*(ctx.gamma + ctx.alpha)*eye(ctx.d) ...
                    + i*ctx.alpha*Gamma{nu};
                Q(rowStart:rowEnd, colStart:colEnd) = block;
            else
                Btilde = computeBtilde(ctx, nu, nuPrime);
                D0_nu_nuPrime = ctx.D{1}(nu, nuPrime);
                Q(rowStart:rowEnd, colStart:colEnd) = ...
                    Btilde + D0_nu_nuPrime * eye(ctx.d);
            end
        end
    end

elseif j == i - 1 && i >= 1  % Subdiagonal block
    for nu = 1:ctx.V
        rowStart = (nu-1)*ctx.d + 1;
        rowEnd = nu*ctx.d;
        colStart = rowStart;
        colEnd = rowEnd;

        block = i*ctx.gamma*eye(ctx.d) + i*ctx.alpha*Bbar{nu};
        Q(rowStart:rowEnd, colStart:colEnd) = block;
    end

elseif j > i && j <= i + ctx.K  % Superdiagonal blocks
    k = j - i;

    for nu = 1:ctx.V
        rowStart = (nu-1)*ctx.d + 1;
        rowEnd = nu*ctx.d;

        for nuPrime = 1:ctx.V
            colStart = (nuPrime-1)*ctx.d + 1;
            colEnd = nuPrime*ctx.d;

            block = zeros(ctx.d, ctx.d);
            for n = 0:ctx.N
                C_nk = computeC(ctx, n, k, nu, nuPrime);
                if ~isempty(C_nk) && any(C_nk(:) ~= 0)
                    nRowStart = getBlockOffset(ctx.T, n);
                    nRowEnd = nRowStart + ctx.T(n+1) - 1;
                    NColStart = getBlockOffset(ctx.T, ctx.N);
                    NColEnd = NColStart + ctx.T(ctx.N+1) - 1;

                    if size(C_nk, 2) == ctx.T(ctx.N+1)
                        block(nRowStart:nRowEnd, NColStart:NColEnd) = C_nk;
                    end
                end
            end
            Q(rowStart:rowEnd, colStart:colEnd) = block;
        end
    end
end
end
