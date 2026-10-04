function result = qsys_dmc(lambda_arr, mu, c, varargin)
% QSYS_DMC Analyzes a D/M/c queue (deterministic interarrivals, exp service).
%
% RESULT = QSYS_DMC(LAMBDA, MU, C) returns time-average performance metrics
% by embedding the system at arrival epochs and integrating over the
% inter-arrival cycle. The death-only sub-generator A[m,m-1] = min(m,c)*MU
% encodes service completions; X_{n+1} = expm(A*s)[X_n+1, :].
%
% Optional name-value parameters:
%   'truncation' - state-space truncation (default auto)
%   'quadSteps'  - trapezoidal-rule steps over a cycle (default 200)
%
% RESULT struct fields:
%   meanQueueLength  - time-average E[N]
%   meanWaitingQueue - time-average Lq = E[(N-c)+]
%   meanWaitingTime  - Wq = Lq / lambda
%   meanSojournTime  - Wq + 1/MU
%   utilization      - rho = lambda / (c*mu)
%   analyzer         - identifier
%
% See also qsys_mapdc, qsys_mdc_crommelin

p = inputParser;
addParameter(p, 'truncation', -1);
addParameter(p, 'quadSteps', 200);
parse(p, varargin{:});
truncation = p.Results.truncation;
quadSteps = p.Results.quadSteps;

if lambda_arr <= 0
    line_error(mfilename, 'Arrival rate must be positive');
end
if mu <= 0
    line_error(mfilename, 'Service rate must be positive');
end
if c < 1
    line_error(mfilename, 'Number of servers must be >= 1');
end
rho = lambda_arr / (c * mu);
if rho >= 1 - 1e-12
    line_error(mfilename, sprintf('Load rho=%g must be strictly less than 1', rho));
end

s = 1.0 / lambda_arr;

if truncation > 0
    nMax = truncation;
else
    nMax = max(200, min(2500, floor(15.0 / (1.0 - rho)) + 200));
end
n = nMax + 1;

% Death-only sub-generator
A = zeros(n);
for m = 0:nMax
    rate = min(m, c) * mu;
    A(m+1, m+1) = -rate;
    if m > 0
        A(m+1, m) = rate;
    end
end

expAs = expm(A * s);
expAdt = expm(A * (s / quadSteps));

% Embedded DTMC at arrival epochs: row X (1-indexed) represents
% "number in system before arrival = X-1". After arrival the state is
% min(X, n-1) which corresponds to row min(X, n-1) + 1 of expAs.
yIdx = min((1:n), n - 1) + 1;
P = zeros(n);
for X = 1:n
    Y = yIdx(X);
    P(X, :) = expAs(Y, :);
end

% Stationary at arrival epochs: solve (P^T - I) pi = 0 with sum pi = 1
M = P' - eye(n);
M(end, :) = 1.0;
b = zeros(n, 1); b(end) = 1.0;
piArr = M \ b;

% Time-average via cycle integration
weightsLq = max((0:nMax) - c, 0);
weightsN = (0:nMax);

LqAtT = zeros(n, quadSteps + 1);
NAtT = zeros(n, quadSteps + 1);
expAt = eye(n);
for k = 0:quadSteps
    if k > 0
        expAt = expAt * expAdt;
    end
    LqAtT(:, k+1) = expAt * weightsLq';
    NAtT(:, k+1) = expAt * weightsN';
end
ts = linspace(0, s, quadSteps + 1);
LqInt = trapz(ts, LqAtT, 2) / s;
NInt = trapz(ts, NAtT, 2) / s;

LqTime = sum(piArr .* LqInt(yIdx));
NTime = sum(piArr .* NInt(yIdx));
Wq = LqTime / lambda_arr;
W = Wq + 1.0 / mu;

result = struct();
result.meanQueueLength = NTime;
result.meanWaitingQueue = LqTime;
result.meanWaitingTime = Wq;
result.meanSojournTime = W;
result.utilization = rho;
result.analyzer = sprintf('Crommelin-DMc:c=%d', c);

end
