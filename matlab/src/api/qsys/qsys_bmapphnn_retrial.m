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
%   'MaxLevel'     - Fixed orbit truncation level. When empty or non-positive
%                    (default) the level is chosen adaptively: it is doubled
%                    until the mass retained at the top level contributes less
%                    than 'TailTolerance' of the mean orbit length. A fixed
%                    level disables the adaptive refinement.
%   'Tolerance'    - Convergence tolerance (default: 1e-10)
%   'TailTolerance'- Relative orbit-truncation error target (default: 1e-6)
%   'MaxDim'       - Cap on the total generator dimension explored by the
%                    adaptive refinement (default: 2e5)
%   'MaxBlockSize' - Cap on the per-level block size V*d (default: 5000).
%                    Exceeding it is an error: the phase-type service order
%                    and the server count make the level block intractable.
%   'Verbose'      - Print progress messages (default: false)
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
%   truncError     - Relative orbit-truncation error estimate at truncLevel
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
addParameter(parser, 'RetrialPolicy', RetrialPolicy.LINEAR, @isnumeric);
addParameter(parser, 'Tolerance', 1e-10, @isnumeric);
addParameter(parser, 'TailTolerance', 1e-6, @isnumeric);
addParameter(parser, 'MaxDim', 2e5, @isnumeric);
addParameter(parser, 'MaxBlockSize', 5000, @isnumeric);
addParameter(parser, 'Verbose', false, @islogical);
parse(parser, varargin{:});

maxLevelParam = parser.Results.MaxLevel;
retrialPolicy = parser.Results.RetrialPolicy;
tol = parser.Results.Tolerance;
tailTol = parser.Results.TailTolerance;
dimMax = parser.Results.MaxDim;
blockMax = parser.Results.MaxBlockSize;
verbose = parser.Results.Verbose;

%% Validate and process inputs

% Reject anything that is not a well-formed Markovian input before any large
% allocation: a malformed (D0,D1) or a NaN service subgenerator otherwise
% propagates silently into the generator and corrupts the solve.
validateRetrialInputs(D, beta, S, N, alpha, gamma, p, R);

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

% Context structure for helper functions
ctx = struct('D', {D}, 'beta', beta, 'S', S, 'S0', S0, 'M', M, 'N', N, ...
    'V', V, 'K', K, 'd', d, 'T', T, 'R', R, 'alpha', alpha, ...
    'gamma', gamma, 'p', p, 'stateMap', {stateMap}, 'retrialPolicy', retrialPolicy);

Vd = V * d;

% A level block of size Vd is dense and is built once per level, so an
% oversized block (high phase-type service order combined with many servers)
% must be rejected rather than attempted.
if Vd > blockMax
    line_error(mfilename, sprintf(['Per-level block size V*d = %d exceeds MaxBlockSize = %d. ' ...
        'The service distribution has %d phases and the station has %d servers, which yields ' ...
        '%d service configurations. Reduce the phase-type order (e.g. fit the service ' ...
        'distribution with fewer phases), reduce the number of servers, or raise ' ...
        '''MaxBlockSize'' if the memory cost is acceptable.'], Vd, round(blockMax), M, N, d));
end

%% Build and solve the system
% see _kb/03-api-layer.md (qsys/ family) for rationale
if ~isempty(maxLevelParam) && maxLevelParam > 0
    truncLevel = round(maxLevelParam);
    pi = solveAtLevel(ctx, truncLevel, Vd, verbose);
    truncError = orbitTruncationError(pi, truncLevel);
else
    truncLevel = max(100, ceil(50 / (1 - min(rho, 0.99))));
    truncError = Inf;
    converged = false;
    while true
        pi = solveAtLevel(ctx, truncLevel, Vd, verbose);
        truncError = orbitTruncationError(pi, truncLevel);
        if truncError <= tailTol
            converged = true;
            break;
        end
        nextLevel = 2 * truncLevel;
        if (nextLevel + 1) * Vd > dimMax
            break;
        end
        if verbose
            fprintf('  Truncation error %.3e > %.3e, refining to level %d\n', ...
                truncError, tailTol, nextLevel);
        end
        truncLevel = nextLevel;
    end
    if ~converged
        line_warning(mfilename, sprintf(['Orbit truncation did not reach the requested accuracy: ' ...
            'residual %.3e > TailTolerance %.3e at level %d (dimension cap MaxDim = %d). ' ...
            'Orbit measures are underestimated; raise ''MaxDim'' or set ''MaxLevel'' explicitly.'], ...
            truncError, tailTol, truncLevel, round(dimMax)));
    end
end

if verbose
    fprintf('  Truncation level: %d (residual %.3e)\n', truncLevel, truncError);
end

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
result.truncError = truncError;
result.analyzer = 'LINE:qsys_bmapphnn_retrial';

if verbose
    fprintf('Solution complete.\n');
end

end

%% ========== Helper Functions ==========

function validateRetrialInputs(D, beta, S, N, alpha, gamma, p, R)
% Reject inputs that are not a well-formed BMAP/PH pair. The generator build
% is driven entirely by these matrices, so a NaN, an Inf, or a non-square
% block would otherwise be written into the generator and only surface as a
% meaningless stationary vector or as an out-of-memory failure.

if ~iscell(D) || isempty(D)
    line_error(mfilename, 'BMAP arrival representation D must be a non-empty cell array {D0,D1,...}.');
end
if numel(D) < 2
    line_error(mfilename, ['BMAP arrival representation D must contain at least {D0,D1}. ' ...
        'A single-element representation is a non-Markovian distribution (e.g. Det or a trace) ' ...
        'and is not admissible in the matrix-analytic retrial engine.']);
end
V = size(D{1}, 1);
for k = 1:numel(D)
    Dk = D{k};
    if ~isnumeric(Dk) || ~ismatrix(Dk) || size(Dk,1) ~= size(Dk,2) || size(Dk,1) ~= V
        line_error(mfilename, sprintf('BMAP matrix D{%d} must be a %dx%d numeric matrix.', k, V, V));
    end
    if any(~isfinite(Dk(:)))
        line_error(mfilename, sprintf(['BMAP matrix D{%d} contains NaN or Inf entries. The arrival ' ...
            'process is disabled or not phase-type representable.'], k));
    end
    if k > 1 && any(Dk(:) < -GlobalConstants.FineTol)
        line_error(mfilename, sprintf('BMAP arrival matrix D{%d} must be non-negative.', k));
    end
end
if any(diag(D{1}) > GlobalConstants.FineTol)
    line_error(mfilename, 'BMAP matrix D0 must have non-positive diagonal entries.');
end
Dsum = zeros(V);
for k = 1:numel(D)
    Dsum = Dsum + D{k};
end
if any(abs(Dsum * ones(V,1)) > sqrt(GlobalConstants.FineTol))
    line_error(mfilename, 'BMAP matrices are inconsistent: sum_k D_k must have zero row sums.');
end

if ~isnumeric(beta) || isempty(beta) || any(~isfinite(beta(:)))
    line_error(mfilename, ['Phase-type service vector beta is empty or contains NaN/Inf. The service ' ...
        'distribution is disabled or not phase-type representable.']);
end
M = size(S, 1);
if ~isnumeric(S) || ~ismatrix(S) || size(S,2) ~= M || numel(beta) ~= M
    line_error(mfilename, sprintf('Phase-type service subgenerator S must be square and conformant with beta (%d phases).', numel(beta)));
end
if any(~isfinite(S(:)))
    line_error(mfilename, ['Phase-type service subgenerator S contains NaN or Inf entries. The service ' ...
        'distribution is disabled or not phase-type representable.']);
end
if any(diag(S) >= 0)
    line_error(mfilename, 'Phase-type service subgenerator S must have strictly negative diagonal entries.');
end
if any(beta(:) < -GlobalConstants.FineTol) || abs(sum(beta(:)) - 1) > sqrt(GlobalConstants.FineTol)
    line_error(mfilename, 'Phase-type service vector beta must be non-negative and sum to one.');
end
if any(-S * ones(M,1) < -GlobalConstants.FineTol)
    line_error(mfilename, 'Phase-type service subgenerator S must have non-negative exit rates.');
end

if ~isscalar(N) || ~isfinite(N) || N < 1 || N ~= round(N)
    line_error(mfilename, 'Number of servers N must be a positive integer.');
end
if ~isscalar(alpha) || ~isfinite(alpha) || alpha < 0
    line_error(mfilename, 'Retrial rate alpha must be a finite non-negative scalar.');
end
if ~isscalar(gamma) || ~isfinite(gamma) || gamma < 0
    line_error(mfilename, 'Orbit impatience rate gamma must be a finite non-negative scalar.');
end
if ~isscalar(p) || ~isfinite(p) || p < 0 || p > 1
    line_error(mfilename, 'Batch rejection probability p must lie in [0,1].');
end
if any(~isfinite(R(:))) || any(R(:) < 0) || any(R(:) > N)
    line_error(mfilename, sprintf('Admission threshold R must lie in [0,%d].', N));
end
end

function err = orbitTruncationError(pi, truncLevel)
% Relative contribution that the truncated tail would add to the mean orbit
% length. Truncation reflects the probability flow that would leave the top
% level back into it, so the mass sitting at the top level bounds the error.
levelMass = sum(pi, 2);
L_orbit = (0:truncLevel) * levelMass;
err = truncLevel * levelMass(end) / max(L_orbit, realmin);
end

function pi = solveAtLevel(ctx, truncLevel, Vd, verbose)
% Build the level-truncated generator and solve pi*Q = 0, pi*e = 1.
%
% The level blocks are level-homogeneous apart from the orbit terms, which
% are linear in the level index: the diagonal block is Qdiag0 + i*Qdiag1
% (retrial and impatience departures from an orbit of size i), the
% subdiagonal block is i*Qsub1 (one of the i orbiting customers succeeds or
% abandons) and the k-th superdiagonal block is level-independent. Building
% those four shapes once and replicating them keeps the assembly linear in
% the truncation level, which the adaptive refinement relies on.
totalDim = (truncLevel + 1) * Vd;

if verbose
    fprintf('Total matrix dimension: %d x %d\n', totalDim, totalDim);
end

% see _kb/03-api-layer.md (qsys/ family) for rationale
ctxGamma = ctx; ctxGamma.alpha = 0;   % impatience only
ctxAlpha = ctx; ctxAlpha.gamma = 0;   % retrials only

Qdiag0 = buildGeneratorLevel(ctx, 0, 0);                       % diagonal block, empty orbit
Qdiag1G = buildGeneratorLevel(ctxGamma, 1, 1) - Qdiag0;        % per-customer impatience increment
Qdiag1A = buildGeneratorLevel(ctxAlpha, 1, 1) - Qdiag0;        % one-unit retrial increment
QsubG = buildGeneratorLevel(ctxGamma, 1, 0);                   % subdiagonal, impatience part
QsubA = buildGeneratorLevel(ctxAlpha, 1, 0);                   % subdiagonal, retrial part
Qsup = cell(1, ctx.K);
for k = 1:ctx.K
    Qsup{k} = buildGeneratorLevel(ctx, 0, k);
end

levels = (0:truncLevel)';

% Retrial weight per level: the orbit size under LINEAR, one whenever the orbit
% is non-empty under CONSTANT.
if ctx.retrialPolicy == RetrialPolicy.CONSTANT
    retrialWeight = double(levels >= 1);
else
    retrialWeight = levels;
end

[rd0, cd0, vd0] = find(Qdiag0);
[rdG, cdG, vdG] = find(Qdiag1G);
[rdA, cdA, vdA] = find(Qdiag1A);
[rsG, csG, vsG] = find(QsubG);
[rsA, csA, vsA] = find(QsubA);

% Diagonal blocks, replicated over all levels
[I, J, X] = replicateBlock(rd0, cd0, vd0, levels, levels, ones(size(levels)), Vd);
% Orbit terms on the diagonal blocks: impatience scales with the orbit size,
% retrials with the policy weight
[I2, J2, X2] = replicateBlock(rdG, cdG, vdG, levels, levels, levels, Vd);
[I2a, J2a, X2a] = replicateBlock(rdA, cdA, vdA, levels, levels, retrialWeight, Vd);
% Subdiagonal blocks (levels 1..truncLevel)
subLevels = levels(levels >= 1);
subWeight = retrialWeight(levels >= 1);
[I3, J3, X3] = replicateBlock(rsG, csG, vsG, subLevels, subLevels - 1, subLevels, Vd);
[I3a, J3a, X3a] = replicateBlock(rsA, csA, vsA, subLevels, subLevels - 1, subWeight, Vd);

I = [I; I2; I2a; I3; I3a];
J = [J; J2; J2a; J3; J3a];
X = [X; X2; X2a; X3; X3a];

for k = 1:ctx.K
    [rk, ck, vk] = find(Qsup{k});
    supLevels = levels(levels <= truncLevel - k);
    if isempty(supLevels) || isempty(rk)
        continue;
    end
    [Ik, Jk, Xk] = replicateBlock(rk, ck, vk, supLevels, supLevels + k, ones(size(supLevels)), Vd);
    I = [I; Ik];
    J = [J; Jk];
    X = [X; Xk];
end

Q = sparse(I, J, X, totalDim, totalDim);

% Ensure rows sum to zero
Q = Q - spdiags(full(sum(Q, 2)), 0, totalDim, totalDim);

% Solve pi * Q = 0, pi * e = 1
if verbose
    fprintf('Solving linear system...\n');
end

Q(:, end) = ones(totalDim, 1);
b = zeros(1, totalDim);
b(end) = 1;

if totalDim > 5000
    if verbose
        fprintf('Using sparse representation\n');
    end
    pi = (Q' \ b')';
else
    pi = b / full(Q);
end

% Reshape to level structure
pi = reshape(pi, Vd, truncLevel + 1)';

% Handle numerical issues
if any(pi(:) < -1e-8)
    line_warning(mfilename, 'Negative probabilities detected, clipping to zero');
end
pi(pi < 0) = 0;
pi = pi / sum(pi(:));  % Renormalize
end

function [I, J, X] = replicateBlock(r, c, v, rowLevels, colLevels, scale, Vd)
% Place a Vd x Vd block pattern (r,c,v) at every (rowLevels, colLevels) pair,
% scaling the entries of the block at position n by scale(n).
if isempty(r) || isempty(rowLevels)
    I = zeros(0,1); J = zeros(0,1); X = zeros(0,1);
    return;
end
r = r(:); c = c(:); v = v(:);
rowLevels = rowLevels(:)'; colLevels = colLevels(:)'; scale = scale(:)';
I = reshape(r + rowLevels * Vd, [], 1);
J = reshape(c + colLevels * Vd, [], 1);
X = reshape(v * scale, [], 1);
end

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
