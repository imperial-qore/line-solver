function result = qsys_phph1(alpha, T, beta, S, varargin)
% QSYS_PHPH1 Analyzes a PH/PH/1 queue using BUTools MMAPPH1FCFS.
%
% RESULT = QSYS_PHPH1(ALPHA, T, BETA, S) analyzes a PH/PH/1 queue with:
%   ALPHA - Arrival PH initial probability vector (1 x n)
%   T     - Arrival PH generator matrix (n x n)
%   BETA  - Service PH initial probability vector (1 x m)
%   S     - Service PH generator matrix (m x m)
%
% The arrival PH is converted to an equivalent MAP representation.
%
% RESULT = QSYS_PHPH1(..., 'numQLMoms', K) computes K queue length moments
% RESULT = QSYS_PHPH1(..., 'numQLProbs', N) computes N queue length probs
% RESULT = QSYS_PHPH1(..., 'numSTMoms', K) computes K sojourn time moments
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
% See also MMAPPH1FCFS, qsys_mapph1, qsys_mapmap1

% Parse optional arguments
p = inputParser;
addParameter(p, 'numQLMoms', 3);
addParameter(p, 'numQLProbs', 100);
addParameter(p, 'numSTMoms', 3);
parse(p, varargin{:});

numQLMoms = p.Results.numQLMoms;
numQLProbs = p.Results.numQLProbs;
numSTMoms = p.Results.numSTMoms;

% Convert arrival PH to MAP representation
% D0 = T (hidden transitions within the PH)
% D1 = (-T)*e * alpha (absorption followed by restart)
[D0, D1] = phToMap(alpha, T);

% Build arrival MMAP structure for BUTools (single class)
D = {D0, D1};

% Service parameters as cell arrays
betaCell = {beta};
SCell = {S};

% Call BUTools solver
[ncMoms, ncDistr, stMoms] = MMAPPH1FCFS(D, betaCell, SCell, ...
    'ncMoms', numQLMoms, 'ncDistr', numQLProbs, 'stMoms', numSTMoms);

% Compute rates
negTinv = inv(-T);
meanInterarrival = alpha * negTinv * ones(size(T, 1), 1);
lambda = 1.0 / meanInterarrival;

negSinv = inv(-S);
meanService = beta * negSinv * ones(size(S, 1), 1);
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

function [D0, D1] = phToMap(alpha, T)
% PHTOMAP Converts a PH distribution to its equivalent MAP representation.
%
% For a PH renewal process, the MAP has:
%   D0 = T (transitions within the PH, no arrival)
%   D1 = t * alpha where t = -T*e (exit rates times restart distribution)

% D0 = T (hidden transitions)
D0 = T;

% t = -T * e (exit rate vector, column)
exitRates = -T * ones(size(T, 1), 1);

% D1 = t * alpha (restart to initial distribution)
D1 = exitRates * alpha;

end
