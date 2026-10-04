function result = qsys_mapg1(D0, D1, serviceMoments, varargin)
% QSYS_MAPG1 Analyzes a MAP/G/1 queue using BUTools MMAPPH1FCFS.
%
% RESULT = QSYS_MAPG1(D0, D1, SERVICEMOMENTS) analyzes a MAP/G/1 queue with:
%   D0             - MAP hidden transition matrix (n x n)
%   D1             - MAP arrival transition matrix (n x n)
%   SERVICEMOMENTS - First k raw moments of service time [E[S], E[S^2], ...]
%                    (k = 2 or 3 for best accuracy)
%
% The general service time is fitted to a Phase-Type distribution using
% moment matching before analysis.
%
% RESULT = QSYS_MAPG1(..., 'numQLMoms', K) computes K queue length moments
% RESULT = QSYS_MAPG1(..., 'numQLProbs', N) computes N queue length probs
% RESULT = QSYS_MAPG1(..., 'numSTMoms', K) computes K sojourn time moments
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
% See also MMAPPH1FCFS, APHFrom3Moments, qsys_mapph1

% Parse optional arguments
p = inputParser;
addParameter(p, 'numQLMoms', 3);
addParameter(p, 'numQLProbs', 100);
addParameter(p, 'numSTMoms', 3);
parse(p, varargin{:});

numQLMoms = p.Results.numQLMoms;
numQLProbs = p.Results.numQLProbs;
numSTMoms = p.Results.numSTMoms;

% Fit service distribution to PH using moment matching
[sigma, S] = fitServiceToPH(serviceMoments);

% Build arrival MMAP structure for BUTools (single class)
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

meanService = serviceMoments(1);
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

function [sigma, S] = fitServiceToPH(moments)
% FITSERVICETOPH Fits a general service time distribution to PH.
%
% Uses moment matching with APHFrom3Moments when 3 moments are available,
% otherwise falls back to simpler approximations.

m1 = moments(1);

if length(moments) >= 3
    % Use APH from 3 moments for best accuracy
    try
        [sigma, S] = APHFrom3Moments(moments(1:3));
    catch
        % Fallback to 2-phase PH
        try
            [sigma, S] = PH2From3Moments(moments(1:3));
        catch
            % Ultimate fallback to exponential
            [sigma, S] = createExponentialPH(m1);
        end
    end
elseif length(moments) >= 2
    % Use 2 moments - create PH(2) with matching mean and variance
    m2 = moments(2);
    cv2 = m2 / (m1 * m1) - 1.0;  % Squared coefficient of variation

    if cv2 <= 0
        % Deterministic or near-deterministic: use Erlang approximation
        k = max(1, round(1.0 / max(cv2, 0.01)));
        [sigma, S] = createErlangPH(m1, k);
    elseif cv2 < 1
        % Hypoexponential: use Erlang-k approximation
        k = max(2, round(1.0 / cv2));
        [sigma, S] = createErlangPH(m1, k);
    elseif cv2 == 1
        % Exponential
        [sigma, S] = createExponentialPH(m1);
    else
        % Hyperexponential: use 2-phase hyperexponential
        [sigma, S] = createHyperexp2PH(m1, cv2);
    end
else
    % Only mean provided - use exponential
    [sigma, S] = createExponentialPH(m1);
end

end

function [sigma, S] = createExponentialPH(mean)
% Creates an exponential PH distribution.
sigma = 1;
S = -1.0 / mean;
end

function [sigma, S] = createErlangPH(mean, k)
% Creates an Erlang-k PH distribution.
mu = k / mean;
sigma = zeros(1, k);
sigma(1) = 1.0;

S = zeros(k, k);
for i = 1:k
    S(i, i) = -mu;
    if i < k
        S(i, i+1) = mu;
    end
end
end

function [sigma, S] = createHyperexp2PH(mean, cv2)
% Creates a 2-phase hyperexponential with matched mean and cv2.

% Balanced means approach
cv = sqrt(cv2);
p = 0.5 * (1.0 + sqrt((cv2 - 1.0) / (cv2 + 1.0)));
lambda1 = 2.0 * p / mean;
lambda2 = 2.0 * (1.0 - p) / mean;

sigma = [p, 1.0 - p];
S = diag([-lambda1, -lambda2]);
end
