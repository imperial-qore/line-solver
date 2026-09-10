function result = qsys_mapph1(D0, D1, sigma, S, varargin)
% QSYS_MAPPH1 Analyzes a MAP/PH/1 queue using BUTools MMAPPH1FCFS.
%
% RESULT = QSYS_MAPPH1(D0, D1, SIGMA, S) analyzes a MAP/PH/1 queue with:
%   D0    - MAP hidden transition matrix (n x n)
%   D1    - MAP arrival transition matrix (n x n)
%   SIGMA - PH service initial probability vector (1 x m)
%   S     - PH service generator matrix (m x m)
%
% RESULT = QSYS_MAPPH1(..., 'numQLMoms', K) computes K queue length moments
% RESULT = QSYS_MAPPH1(..., 'numQLProbs', N) computes N queue length probs
% RESULT = QSYS_MAPPH1(..., 'numSTMoms', K) computes K sojourn time moments
%
% Returns a struct with fields:
%   meanQueueLength   - Mean number of customers in system
%   meanWaitingTime   - Mean waiting time in queue
%   meanSojournTime   - Mean sojourn time (waiting + service)
%   utilization       - Server utilization
%   queueLengthDist   - Queue length distribution P(Q=n)
%   queueLengthMoments- Raw moments of queue length
%   sojournTimeMoments- Raw moments of sojourn time
%   analyzer          - Name of analyzer used
%
% See also MMAPPH1FCFS, qsys_mapmap1, qsys_phph1

% Parse optional arguments
p = inputParser;
addParameter(p, 'numQLMoms', 3);
addParameter(p, 'numQLProbs', 100);
addParameter(p, 'numSTMoms', 3);
parse(p, varargin{:});

numQLMoms = p.Results.numQLMoms;
numQLProbs = p.Results.numQLProbs;
numSTMoms = p.Results.numSTMoms;

% Build MMAP structure for BUTools (single class)
D = {D0, D1};

% Service parameters as cell arrays
sigmaCell = {sigma};
SCell = {S};

% Call BUTools solver
[ncMoms, ncDistr, stMoms] = MMAPPH1FCFS(D, sigmaCell, SCell, ...
    'ncMoms', numQLMoms, 'ncDistr', numQLProbs, 'stMoms', numSTMoms);

% Compute utilization from arrival and service rates
theta = ctmc_solve(D0 + D1);
lambda = sum(theta * D1);

negSinv = inv(-S);
meanService = sigma * negSinv * ones(size(S,1), 1);
mu = 1.0 / meanService;
rho = lambda / mu;

% Extract results
if ~isempty(ncMoms)
    meanQL = ncMoms(1);
else
    meanQL = 0;
end

if ~isempty(stMoms)
    meanST = stMoms(1);
else
    meanST = 0;
end

% Waiting time = sojourn time - service time
meanWT = max(0, meanST - meanService);

% Build result struct
result = struct();
result.meanQueueLength = meanQL;
result.meanWaitingTime = meanWT;
result.meanSojournTime = meanST;
result.utilization = rho;
result.queueLengthDist = ncDistr;
result.queueLengthMoments = ncMoms;
result.sojournTimeMoments = stMoms;
result.analyzer = 'BUTools:MMAPPH1FCFS';

end
