function result = qsys_mapmc(D0, D1, mu, c, varargin)
% QSYS_MAPMC Analyzes a MAP/M/c queue using Q-MAM.
%
% RESULT = QSYS_MAPMC(D0, D1, MU, C) analyzes a MAP/M/c queue with:
%   D0 - MAP hidden transition matrix (n x n)
%   D1 - MAP arrival transition matrix (n x n)
%   MU - Exponential service rate
%   C  - Number of servers
%
% RESULT = QSYS_MAPMC(..., 'maxNumComp', N) sets max queue length probs (default 500)
%
% Returns a struct with fields:
%   meanQueueLength   - Mean number of customers in system
%   meanWaitingTime   - Mean waiting time in queue
%   meanSojournTime   - Mean sojourn time (waiting + service)
%   utilization       - Server utilization (per server)
%   queueLengthDist   - Queue length distribution P(Q=n)
%   waitingTimePH     - Struct with alpha and T for waiting time PH
%   analyzer          - Name of analyzer used
%
% See also Q_CT_MAP_M_C, qsys_mapph1, qsys_mapm1

% Parse optional arguments
p = inputParser;
addParameter(p, 'maxNumComp', 500);
parse(p, varargin{:});

maxNumComp = p.Results.maxNumComp;

% Call Q-MAM solver
[ql, wait_alpha, Smat] = Q_CT_MAP_M_C(D0, D1, mu, c, 'MaxNumComp', maxNumComp);

% Compute arrival rate from MAP
theta = ctmc_solve(D0 + D1);
lambda = sum(theta * D1);

% Compute utilization
rho = lambda / (mu * c);

% Compute mean queue length from distribution
meanQL = 0;
for i = 1:length(ql)
    meanQL = meanQL + (i - 1) * ql(i);
end

% Compute mean waiting time from PH representation
if ~isempty(wait_alpha) && ~isempty(Smat)
    % Mean of PH distribution = alpha * (-T)^{-1} * e
    nonzeroIdx = find(wait_alpha > 0);
    if ~isempty(nonzeroIdx)
        wait_alpha_nz = wait_alpha(nonzeroIdx);
        negSinv = inv(-Smat);
        meanWT = wait_alpha_nz * negSinv * ones(size(Smat, 1), 1);
    else
        meanWT = 0;
    end
else
    meanWT = 0;
end

% Mean service time
meanService = 1.0 / mu;

% Mean sojourn time = waiting + service
meanST = meanWT + meanService;

% Build result struct
result = struct();
result.meanQueueLength = meanQL;
result.meanWaitingTime = meanWT;
result.meanSojournTime = meanST;
result.utilization = rho;
result.queueLengthDist = ql;
result.waitingTimePH = struct('alpha', wait_alpha, 'T', Smat);
result.analyzer = sprintf('Q-MAM:MAP/M/%d', c);

end
