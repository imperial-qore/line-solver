function result = qsys_mapdc(D0, D1, s, c, varargin)
% QSYS_MAPDC Analyzes a MAP/D/c queue using Q-MAM.
%
% RESULT = QSYS_MAPDC(D0, D1, S, C) analyzes a MAP/D/c queue with:
%   D0 - MAP hidden transition matrix (n x n)
%   D1 - MAP arrival transition matrix (n x n)
%   S  - Deterministic service time (positive scalar)
%   C  - Number of servers (positive integer)
%
% RESULT = QSYS_MAPDC(..., 'maxNumComp', N) sets max queue length components (default 1000)
% RESULT = QSYS_MAPDC(..., 'numSteps', K) sets waiting time distribution granularity (default 1)
%
% Returns a struct with fields:
%   meanQueueLength   - Mean number of customers in system
%   meanWaitingTime   - Mean waiting time in queue
%   meanSojournTime   - Mean sojourn time (waiting + service)
%   utilization       - Server utilization (per server)
%   queueLengthDist   - Queue length distribution P(Q=n)
%   waitingTimeDist   - Waiting time CDF at discrete points
%   analyzer          - Name of analyzer used
%
% The waiting time distribution is evaluated at points {0, s/numSteps, 2*s/numSteps, ...}
% where s is the deterministic service time and numSteps controls the granularity.
%
% See also Q_CT_MAP_D_C, qsys_mapmc, qsys_mapd1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Parse optional arguments
p = inputParser;
addParameter(p, 'maxNumComp', 1000);
addParameter(p, 'numSteps', 1);
parse(p, varargin{:});

maxNumComp = p.Results.maxNumComp;
numSteps = p.Results.numSteps;

% Input validation
if ~isnumeric(s) || ~isreal(s) || s <= 0
    line_error(mfilename, 'Service time s must be a positive real number');
end
if ~isnumeric(c) || ~isreal(c) || c < 1 || floor(c) ~= c
    line_error(mfilename, 'Number of servers c must be a positive integer');
end

% Call Q-MAM solver
[ql, w] = Q_CT_MAP_D_C(D0, D1, s, c, 'MaxNumComp', maxNumComp, 'NumSteps', numSteps);

% Compute arrival rate from MAP
theta = ctmc_solve(D0 + D1);
lambda = sum(theta * D1);

% Compute utilization (note: s is service time, not rate)
rho = lambda * s / c;

% Compute mean queue length from distribution
% ql(i) = Prob[(i-1) customers in the queue]
n = length(ql);
meanQL = sum((0:n-1) .* ql(:)');

% Compute mean waiting time from CDF
% w(i) = Prob[waiting time <= (i-1)*s/numSteps]
% Mean = integral of survival function = sum((1 - CDF) * stepSize)
if ~isempty(w) && length(w) > 1
    stepSize = s / numSteps;
    % Survival function at each point (except the last)
    survivalFn = 1 - w(1:end-1);
    meanWT = sum(survivalFn(:)) * stepSize;
else
    % Fallback: estimate from queue length using Little's law approximation
    meanWT = max(0, (meanQL - rho * c) / lambda);
end

% Mean service time (deterministic)
meanService = s;

% Mean sojourn time = waiting + service
meanST = meanWT + meanService;

% Build result struct
result = struct();
result.meanQueueLength = meanQL;
result.meanWaitingTime = meanWT;
result.meanSojournTime = meanST;
result.utilization = rho;
result.queueLengthDist = ql;
result.waitingTimeDist = w;
result.analyzer = sprintf('Q-MAM:MAP/D/%d', c);

end
