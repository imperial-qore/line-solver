function result = qsys_bmapm1(D, mu, varargin)
% QSYS_BMAPM1 Analyzes a BMAP/M/1 queue by the matrix-analytic (M/G/1-type) method.
%
% RESULT = QSYS_BMAPM1(D, MU) analyzes a single-server queue fed by a batch
% Markovian arrival process and with exponential service of rate MU.
%
% Inputs:
%   D  - cell array {D0, D1, ..., DK} of BMAP matrices. D0 carries the hidden
%        transitions, Dk (k >= 1) the transitions that release a batch of k
%        customers.
%   MU - exponential service rate.
%
% Optional parameters:
%   'Uniformization' - uniformization constant q used to randomize the
%                      generator into a discrete-time M/G/1-type chain. It must
%                      dominate every total outflow rate; by default it is
%                      chosen as max_i(-D0(i,i)) + mu, rounded up.
%   'MaxIter'        - maximum functional iterations for G (default 10000)
%   'Tolerance'      - convergence tolerance for G (default 1e-12)
%   'MaxLevel'       - level truncation used for the queue-length distribution
%                      (default: adaptive, see qsys_bmapphnn_retrial)
%   'TailTolerance'  - relative truncation target for the level distribution
%                      (default 1e-10)
%
% Beyond the usual performance measures the result exposes the intermediate
% matrix-analytic quantities themselves, so that the algorithm can be inspected
% and taught rather than only its output:
%
%   theta        - stationary vector of the BMAP phase process, sum_k D_k
%   lambda       - mean arrival rate, theta * sum_k k*D_k * e
%   rho          - offered load lambda/mu
%   q            - uniformization constant actually used
%   A0, A1, Bk   - randomized blocks: A0 = (mu/q)I is a service completion
%                  (level down by one), A1 = (1/q)(D0 - mu*I) + I keeps the
%                  level, Bk{k} = (1/q)D_k raises it by k
%   B0           - boundary local block (1/q)D0 + I, used at level 0 where no
%                  service can complete
%   A            - A0 + A1 + sum_k Bk{k}, the phase process of the chain
%   alpha        - stationary vector of A
%   G            - minimal non-negative solution of
%                  G = A0 + A1*G + sum_k Bk{k}*G^(k+1)
%   drift        - alpha*(sum_k k*Bk{k})*e - alpha*A0*e. The queue is stable
%                  iff this is strictly negative
%   decayRate    - geometric decay rate of the level probabilities, measured as
%                  the limiting ratio pi_(n+1)/pi_n. Reported rather than
%                  derived from a spectral convention so that it is unambiguous
%   levelProb    - level probabilities pi_n as rows (level 0 first)
%   pi0          - probability the system is empty (equals 1-rho exactly)
%
% Example:
%   % Example 6.4 of Bolch et al.
%   D0 = [-2, 1/2; 1/3, -3]; D1 = [1/4, 1/2; 1/3, 1]; D2 = [1/4, 1/2; 1, 1/3];
%   result = qsys_bmapm1({D0, D1, D2}, 11);
%
% See also qsys_mapm1, qsys_mapph1, qsys_bmapphnn_retrial
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

parser = inputParser;
addParameter(parser, 'Uniformization', [], @isnumeric);
addParameter(parser, 'MaxIter', 10000, @isnumeric);
addParameter(parser, 'Tolerance', 1e-12, @isnumeric);
addParameter(parser, 'MaxLevel', [], @isnumeric);
addParameter(parser, 'TailTolerance', 1e-10, @isnumeric);
parse(parser, varargin{:});

qParam = parser.Results.Uniformization;
maxIter = parser.Results.MaxIter;
tol = parser.Results.Tolerance;
maxLevelParam = parser.Results.MaxLevel;
tailTol = parser.Results.TailTolerance;

%% Validate inputs
if ~iscell(D) || numel(D) < 2
    line_error(mfilename, 'The BMAP must be given as a cell array {D0,D1,...,DK} with at least D0 and D1.');
end
V = size(D{1}, 1);
for k = 1:numel(D)
    Dk = D{k};
    if ~isnumeric(Dk) || size(Dk,1) ~= V || size(Dk,2) ~= V || any(~isfinite(Dk(:)))
        line_error(mfilename, sprintf('BMAP matrix D{%d} must be a finite %dx%d matrix.', k, V, V));
    end
    if k > 1 && any(Dk(:) < -GlobalConstants.FineTol)
        line_error(mfilename, sprintf('BMAP arrival matrix D{%d} must be non-negative.', k));
    end
end
Dsum = zeros(V);
for k = 1:numel(D)
    Dsum = Dsum + D{k};
end
if any(abs(Dsum * ones(V,1)) > sqrt(GlobalConstants.FineTol))
    line_error(mfilename, 'BMAP matrices are inconsistent: sum_k D_k must have zero row sums.');
end
if ~isscalar(mu) || ~isfinite(mu) || mu <= 0
    line_error(mfilename, 'The service rate mu must be a finite positive scalar.');
end

K = numel(D) - 1;

%% Arrival characterization
theta = ctmc_solve(Dsum);
theta = theta(:)';
sumKDk = zeros(V);
for k = 1:K
    sumKDk = sumKDk + k * D{k+1};
end
lambda = theta * sumKDk * ones(V,1);
rho = lambda / mu;

%% Randomization (uniformization) into a discrete-time M/G/1-type chain
if isempty(qParam)
    q = max(-diag(D{1})) + mu;
else
    q = qParam;
end
if q < max(-diag(D{1})) + mu - GlobalConstants.FineTol
    line_error(mfilename, sprintf(['The uniformization constant q = %g does not dominate the total outflow rate ' ...
        '%g; the randomized chain would have negative entries.'], q, max(-diag(D{1})) + mu));
end

A0 = (mu/q) * eye(V);                        % level down by one: service completion
A1 = (1/q) * (D{1} - mu*eye(V)) + eye(V);    % level unchanged
B0 = (1/q) * D{1} + eye(V);                  % level 0: no service can complete
Bk = cell(1, K);
for k = 1:K
    Bk{k} = (1/q) * D{k+1};                  % level up by k
end

A = A0 + A1;
for k = 1:K
    A = A + Bk{k};
end
alpha = dtmc_solve(A);
alpha = alpha(:)';

%% Matrix G: minimal non-negative solution of the M/G/1-type equation
G = zeros(V);
converged = false;
for iter = 1:maxIter
    Gpow = G;
    Gnew = A0 + A1*G;
    for k = 1:K
        Gpow = Gpow * G;          % G^(k+1)
        Gnew = Gnew + Bk{k} * Gpow;
    end
    if max(abs(Gnew(:) - G(:))) < tol
        G = Gnew;
        converged = true;
        break;
    end
    G = Gnew;
end
if ~converged
    line_warning(mfilename, sprintf(['The functional iteration for G did not converge to %g in %d iterations ' ...
        '(last change %g). The queue may be unstable.'], tol, maxIter, max(abs(Gnew(:) - G(:)))));
end

%% Stability drift
upDrift = zeros(V);
for k = 1:K
    upDrift = upDrift + k * Bk{k};
end
drift = alpha * upDrift * ones(V,1) - alpha * A0 * ones(V,1);

%% Level probabilities of the continuous-time chain
% Built from the level-truncated generator: the level blocks are homogeneous
% above the boundary, so a single truncated solve gives the whole distribution
% up to a residual that is refined until negligible.
if ~isempty(maxLevelParam) && maxLevelParam > 0
    levelMax = round(maxLevelParam);
    levelProb = solveLevels(D, mu, V, K, levelMax);
    truncError = levelTailError(levelProb, levelMax);
else
    levelMax = max(50, ceil(20 / max(1 - min(rho, 0.999), eps)));
    truncError = Inf;
    while true
        levelProb = solveLevels(D, mu, V, K, levelMax);
        truncError = levelTailError(levelProb, levelMax);
        if truncError <= tailTol || (2*levelMax + 1) * V > 2e5
            break;
        end
        levelMax = 2 * levelMax;
    end
    if truncError > tailTol
        line_warning(mfilename, sprintf(['The level distribution did not reach the requested accuracy: residual ' ...
            '%.3e > TailTolerance %.3e at level %d.'], truncError, tailTol, levelMax));
    end
end

levelMass = sum(levelProb, 2);
% Measured decay rate: the ratio settles geometrically, so read it where the
% mass is still numerically meaningful rather than at the truncation boundary.
usable = find(levelMass > 1e-12, 1, 'last');
if isempty(usable) || usable < 3
    decayRate = NaN;
else
    ref = max(2, floor(usable/2));
    decayRate = levelMass(ref+1) / levelMass(ref);
end

meanQueueLength = (0:(size(levelProb,1)-1)) * levelMass;

%% Result
result = struct();
result.theta = theta;
result.lambda = lambda;
result.rho = rho;
result.q = q;
result.A0 = A0;
result.A1 = A1;
result.B0 = B0;
result.Bk = {Bk};
result.A = A;
result.alpha = alpha;
result.G = G;
result.drift = drift;
result.decayRate = decayRate;
result.levelProb = levelProb;
result.pi0 = levelMass(1);
result.meanQueueLength = meanQueueLength;
result.utilization = rho;
result.throughput = lambda;
result.truncLevel = size(levelProb,1) - 1;
result.truncError = truncError;
result.analyzer = 'LINE:qsys_bmapm1';
end

function err = levelTailError(levelProb, levelMax)
% Relative contribution the truncated tail would add to the mean level.
levelMass = sum(levelProb, 2);
meanLevel = (0:levelMax) * levelMass;
err = levelMax * levelMass(end) / max(meanLevel, realmin);
end

function levelProb = solveLevels(D, mu, V, K, levelMax)
% Level-truncated CTMC generator of the BMAP/M/1 queue and its stationary
% distribution. Level n holds n customers in the system; the phase is the BMAP
% state. Service fires only above level 0.
totalDim = (levelMax + 1) * V;
levels = (0:levelMax)';

[r0, c0, v0] = find(D{1});
[rs, cs, vs] = find(mu * eye(V));

% Local blocks: D0 on every level, minus the service rate above level 0 (the
% service outflow is put back on the diagonal by the row-sum correction).
[I, J, X] = tileBlock(r0, c0, v0, levels, levels, V);
% Service: level n -> n-1 for n >= 1
sub = levels(levels >= 1);
[Is, Js, Xs] = tileBlock(rs, cs, vs, sub, sub - 1, V);
I = [I; Is]; J = [J; Js]; X = [X; Xs];
% Batch arrivals: level n -> n+k
for k = 1:K
    [rk, ck, vk] = find(D{k+1});
    if isempty(rk)
        continue
    end
    up = levels(levels <= levelMax - k);
    if isempty(up)
        continue
    end
    [Ik, Jk, Xk] = tileBlock(rk, ck, vk, up, up + k, V);
    I = [I; Ik]; J = [J; Jk]; X = [X; Xk];
end

Q = sparse(I, J, X, totalDim, totalDim);
Q = Q - spdiags(full(sum(Q, 2)), 0, totalDim, totalDim);

Q(:, end) = ones(totalDim, 1);
b = zeros(1, totalDim);
b(end) = 1;
if totalDim > 5000
    pi = (Q' \ b')';
else
    pi = b / full(Q);
end

levelProb = reshape(pi, V, levelMax + 1)';
levelProb(levelProb < 0) = 0;
levelProb = levelProb / sum(levelProb(:));
end

function [I, J, X] = tileBlock(r, c, v, rowLevels, colLevels, V)
% Place a V x V block pattern at every (rowLevels, colLevels) pair.
if isempty(r) || isempty(rowLevels)
    I = zeros(0,1); J = zeros(0,1); X = zeros(0,1);
    return
end
r = r(:); c = c(:); v = v(:);
rowLevels = rowLevels(:)'; colLevels = colLevels(:)';
I = reshape(r + rowLevels * V, [], 1);
J = reshape(c + colLevels * V, [], 1);
X = repmat(v, numel(rowLevels), 1);
end
