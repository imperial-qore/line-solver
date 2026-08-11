function result = qsys_mapmap1(C0, C1, D0, D1, varargin)
% QSYS_MAPMAP1 Analyzes a MAP/MAP/1 queue using BUTools MMAPPH1FCFS.
%
% RESULT = QSYS_MAPMAP1(C0, C1, D0, D1) analyzes a MAP/MAP/1 queue with:
%   C0    - Arrival MAP hidden transition matrix (n x n)
%   C1    - Arrival MAP arrival transition matrix (n x n)
%   D0    - Service MAP hidden transition matrix (m x m)
%   D1    - Service MAP observable transition matrix (m x m)
%
% The service MAP is converted to an equivalent PH representation.
%
% RESULT = QSYS_MAPMAP1(..., 'numQLMoms', K) computes K queue length moments
% RESULT = QSYS_MAPMAP1(..., 'numQLProbs', N) computes N queue length probs
% RESULT = QSYS_MAPMAP1(..., 'numSTMoms', K) computes K sojourn time moments
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
% See also MMAPPH1FCFS, qsys_mapph1, qsys_phph1

% Parse optional arguments
p = inputParser;
addParameter(p, 'numQLMoms', 3);
addParameter(p, 'numQLProbs', 100);
addParameter(p, 'numSTMoms', 3);
parse(p, varargin{:});

numQLMoms = p.Results.numQLMoms;
numQLProbs = p.Results.numQLProbs;
numSTMoms = p.Results.numSTMoms;

% Build arrival MMAP structure for BUTools (single class)
D = {C0, C1};

% Convert service MAP to PH representation
[sigma, S] = mapToPh(D0, D1);

% Service parameters as cell arrays
sigmaCell = {sigma};
SCell = {S};

% Call BUTools solver
[ncMoms, ncDistr, stMoms] = MMAPPH1FCFS(D, sigmaCell, SCell, ...
    'ncMoms', numQLMoms, 'ncDistr', numQLProbs, 'stMoms', numSTMoms);

% Compute utilization from arrival and service rates
thetaArr = ctmc_solve(C0 + C1);
lambda = sum(thetaArr * C1);

thetaSvc = ctmc_solve(D0 + D1);
mu = sum(thetaSvc * D1);
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

% Mean service time from PH
negSinv = inv(-S);
meanService = sigma * negSinv * ones(size(S,1), 1);

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

function [sigma, S] = mapToPh(D0, D1)
% MAPTOPH Converts a MAP to its equivalent PH representation.
%
% The PH representation captures the service time distribution
% of the MAP, with initial distribution derived from the
% stationary distribution weighted by observable transition rates.

% Stationary distribution of the MAP
theta = ctmc_solve(D0 + D1);

% Initial distribution for PH is proportional to
% how customers enter the service process
sigma = theta * D1;
total = sum(sigma);
if total > 0
    sigma = sigma / total;
else
    % Fallback to uniform if no arrivals
    sigma = ones(1, size(D0, 1)) / size(D0, 1);
end

% The PH generator is the hidden transition matrix
S = D0;

end
