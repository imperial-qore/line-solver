function result = qsys_mmck(lambda, mu, c, K)
% QSYS_MMCK Exact closed-form analysis of an M/M/c/K queue.
%
% RESULT = QSYS_MMCK(LAMBDA, MU, C, K) analyzes an M/M/c/K queue with:
%   LAMBDA - Poisson arrival rate (positive scalar)
%   MU     - Exponential service rate per server (positive scalar)
%   C      - Number of servers (positive integer)
%   K      - System capacity, total jobs allowed (integer, K >= C)
%
% Stationary distribution (Erlang-B/C truncated form):
%   a = lambda/mu, rho = a/c
%   p_n = (a^n / n!) * p_0                    for 0 <= n <= c
%   p_n = (a^c / c!) * rho^(n-c) * p_0        for c <= n <= K
%   p_0 = 1 / [ sum_{n=0}^{c-1} a^n/n! + (a^c/c!) * sum_{n=0}^{K-c} rho^n ]
%
% Returns a struct with fields:
%   meanQueueLength    - Mean number of jobs in system, L
%   meanQueueLengthQ   - Mean number waiting in queue, Lq
%   meanWaitingTime    - Mean waiting time in queue, Wq (Little)
%   meanSojournTime    - Mean sojourn time, W = Wq + 1/mu
%   utilization        - Per-server utilization, lambda_eff/(c*mu)
%   throughput         - Effective throughput, lambda*(1 - p_K)
%   lossProbability    - Blocking probability, p_K
%   queueLengthDist    - Distribution as 1x(K+1) row: [p_0, p_1, ..., p_K]
%   analyzer           - Identifier string
%
% Examples:
%   r = qsys_mmck(1, 1, 2, 4);   % M/M/2/4 with rho = 0.5
%   r.meanQueueLength            % 1.1304
%   r.lossProbability            % 0.0435
%
% See also QSYS_MM1K_LOSS, QSYS_MMK, QSYS_MAPDC
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Input validation
if ~isnumeric(lambda) || ~isscalar(lambda) || ~isreal(lambda) || lambda <= 0
    line_error(mfilename, 'lambda must be a positive real scalar');
end
if ~isnumeric(mu) || ~isscalar(mu) || ~isreal(mu) || mu <= 0
    line_error(mfilename, 'mu must be a positive real scalar');
end
if ~isnumeric(c) || ~isscalar(c) || ~isreal(c) || c < 1 || floor(c) ~= c
    line_error(mfilename, 'c must be a positive integer');
end
if ~isnumeric(K) || ~isscalar(K) || ~isreal(K) || K < c || floor(K) ~= K
    line_error(mfilename, 'K must be an integer >= c');
end

a = lambda / mu;
rho = a / c;

% Build unnormalized stationary distribution
p = zeros(1, K+1);
% Levels 0..c-1: Erlang-B partial sum terms
for n = 0:(c-1)
    p(n+1) = a^n / factorial(n);
end
% Levels c..K: geometric tail in rho beyond level c
ac_over_cfact = a^c / factorial(c);
for n = c:K
    p(n+1) = ac_over_cfact * rho^(n-c);
end
% Normalize
S = sum(p);
if ~isfinite(S) || S <= 0
    line_error(mfilename, 'Stationary distribution failed to normalize (numerical overflow?)');
end
p = p / S;

% Performance metrics
levels = 0:K;
L = levels * p(:);                    % Mean # in system
p_K = p(K+1);                         % Blocking probability
lambdaEff = lambda * (1 - p_K);       % Effective throughput
n_waiting = max(0, levels - c);       % Jobs waiting at each level
Lq = n_waiting * p(:);                % Mean # in queue
util = lambdaEff / (c * mu);          % Per-server utilization
if lambdaEff > 0
    Wq = Lq / lambdaEff;
    W = L / lambdaEff;
else
    Wq = 0;
    W = 0;
end

result = struct();
result.meanQueueLength  = L;
result.meanQueueLengthQ = Lq;
result.meanWaitingTime  = Wq;
result.meanSojournTime  = W;
result.utilization      = util;
result.throughput       = lambdaEff;
result.lossProbability  = p_K;
result.queueLengthDist  = p;
result.analyzer         = sprintf('exact:M/M/%d/%d', c, K);
end
