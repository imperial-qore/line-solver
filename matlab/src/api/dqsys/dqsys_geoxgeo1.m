function result = dqsys_geoxgeo1(a, beta, s, convention)
% QSYS_GEOXGEO1 Exact analysis of a discrete-time Geo^X/Geo/1 batch queue.
%
% RESULT = QSYS_GEOXGEO1(A, BETA, S) analyzes a slotted single-server queue
% with batch arrivals:
%   A    - per-slot probability that a batch arrives (0 < A <= 1)
%   BETA - batch-size geometric parameter; the batch is supported on
%          {1,2,...} with mean 1/BETA (0 < BETA <= 1)
%   S    - per-slot service completion probability (0 < S <= 1)
% Stability requires LAMBDA = A/BETA < S.
%
% RESULT = QSYS_GEOXGEO1(A, BETA, S, CONVENTION) selects the observation
% epoch, 'LAS_DA' (default) or 'EAS'; see QSYS_GEOGEO1 for what they mean.
%
% With A(z) = 1-A+A*X(z) the pgf of the number of jobs arriving in one slot,
% the slot-boundary content obeys X(t+1) = X(t) - D(t) + A(t) with the
% departure resolved first, giving
%   P(z) = p_0 S (z-1) A(z) / ( z - A(z)(S+(1-S)z) ),   p_0 = 1 - LAMBDA/S
% and, on differentiating at z = 1,
%   E[N] = LAMBDA + ( A E[X(X-1)]/2 + LAMBDA(1-S) ) / (S - LAMBDA)
% The batch enters only through its first two factorial moments. For the
% geometric batch, E[X] = 1/BETA and E[X(X-1)] = 2(1-BETA)/BETA^2. At BETA = 1
% the batch is always a single job and the result equals QSYS_GEOGEO1.
%
% Returns a struct with fields:
%   convention        - Observation epoch used, 'LAS_DA' or 'EAS'
%   batchArrivalProb  - Per-slot probability that a batch arrives, A
%   batchMean         - Mean batch size E[X]
%   batchSecondFactorialMoment - E[X(X-1)]
%   serviceProb       - Per-slot service completion probability S
%   arrivalRate       - Jobs arriving per slot, LAMBDA = A E[X]
%   throughput        - Departures per slot, LAMBDA
%   utilization       - Fraction of slots with the server serving, LAMBDA/S
%   boundaryEmptyProb - Empty probability AT THE SLOT BOUNDARY, 1 - LAMBDA/S
%   meanQueueLength   - Mean number of jobs in the system
%   meanWaitingQueue  - Mean number of jobs waiting, i.e. not in service
%   meanSojournTime   - Mean sojourn time in slots, per job (not per batch)
%   meanWaitingTime   - Mean waiting time in slots, per job
%   meanServiceTime   - Mean service time, 1/S under LAS_DA, (1-S)/S under EAS
%   pgf               - Function handle P(z, Az) taking the slot-arrival pgf
%                       value A(z) at the same z; defined for 0 < z <= 1
%   analyzer          - Identifier string
%
% No pmf is returned: for a general batch law the stationary distribution has
% no elementary closed form, so only the generating function is exact.
% Likewise BOUNDARYEMPTYPROB is deliberately epoch independent; the empty
% probability at the EAS epoch is p_0 + S p_1 and has no elementary form.
%
% Examples:
%   r = dqsys_geoxgeo1(0.1, 0.5, 0.9);   % lambda = 0.2
%   r.arrivalRate                        % 0.2
%   r.meanQueueLength                    % 0.5143
%
% Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
% Springer 2001, chapter 6 for the discrete-time batch framework; the
% single-node pgf above is the standard discrete-time M/G/1-type derivation.
%
% See also QSYS_GEOGEO1, QSYS_MXM1
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(convention)
    convention = 'LAS_DA';
end
if ~isnumeric(beta) || ~isscalar(beta) || ~isreal(beta) || beta <= 0 || beta > 1
    line_error(mfilename, 'beta must be a real scalar in (0,1]');
end
batchMean = 1 / beta;
batchSecondFactorial = 2 * (1 - beta) / (beta ^ 2);
result = dqsys_geoxgeo1_moments(a, batchMean, batchSecondFactorial, s, convention);
end

function result = dqsys_geoxgeo1_moments(a, batchMean, batchSecondFactorial, s, convention)
% Geo^X/Geo/1 for an arbitrary batch law given by its first two factorial
% moments. Kept local because the geometric form above is the documented
% entry point; call QSYS_GEOXGEO1 with the matching beta for other laws, or
% lift this to its own file if a non-geometric batch is needed directly.

if ~ischar(convention) && ~isstring(convention)
    line_error('dqsys_geoxgeo1', 'convention must be the string ''LAS_DA'' or ''EAS''');
end
convention = upper(char(convention));
if ~strcmp(convention, 'LAS_DA') && ~strcmp(convention, 'EAS')
    line_error('dqsys_geoxgeo1', 'convention must be ''LAS_DA'' or ''EAS''');
end
if ~isnumeric(a) || ~isscalar(a) || ~isreal(a) || a <= 0 || a > 1
    line_error('dqsys_geoxgeo1', 'a must be a real scalar in (0,1]');
end
if ~isnumeric(s) || ~isscalar(s) || ~isreal(s) || s <= 0 || s > 1
    line_error('dqsys_geoxgeo1', 's must be a real scalar in (0,1]');
end
if batchMean < 1
    line_error('dqsys_geoxgeo1', 'mean batch size must be at least 1: a batch that arrives carries at least one job');
end
if batchSecondFactorial < 0
    line_error('dqsys_geoxgeo1', 'E[X(X-1)] must be non-negative');
end
% E[X^2] = E[X(X-1)] + E[X] >= E[X]^2, so a smaller second factorial moment
% describes no random variable at all.
minSecondFactorial = batchMean ^ 2 - batchMean;
if batchSecondFactorial < minSecondFactorial - 1e-9 * max(1, minSecondFactorial)
    line_error('dqsys_geoxgeo1', 'E[X(X-1)] is below E[X]^2-E[X], so the batch moments describe no random variable');
end

lambda = a * batchMean;
if lambda >= s
    line_error('dqsys_geoxgeo1', 'load lambda/s must be strictly less than 1');
end

rho = lambda / s;
boundaryEmptyProb = 1 - rho;

% P'(1) from the generating function.
meanAtBoundary = lambda + (a * batchSecondFactorial / 2 + lambda * (1 - s)) / (s - lambda);

meanSojournAtBoundary = meanAtBoundary / lambda;
meanWaitingTime = meanSojournAtBoundary - 1 / s;
meanWaitingQueue = lambda * meanWaitingTime;

if strcmp(convention, 'LAS_DA')
    meanQueueLength = meanAtBoundary;
    meanSojournTime = meanSojournAtBoundary;
    meanServiceTime = 1 / s;
else
    % One departure earlier: the epoch drops exactly the departures of the
    % slot, whose rate is lambda, hence one slot of sojourn.
    meanQueueLength = meanAtBoundary - lambda;
    meanSojournTime = meanSojournAtBoundary - 1;
    meanServiceTime = (1 - s) / s;
end

result = struct();
result.convention = convention;
result.batchArrivalProb = a;
result.batchMean = batchMean;
result.batchSecondFactorialMoment = batchSecondFactorial;
result.serviceProb = s;
result.arrivalRate = lambda;
result.throughput = lambda;
result.utilization = rho;
result.boundaryEmptyProb = boundaryEmptyProb;
result.meanQueueLength = meanQueueLength;
result.meanWaitingQueue = meanWaitingQueue;
result.meanSojournTime = meanSojournTime;
result.meanWaitingTime = meanWaitingTime;
result.meanServiceTime = meanServiceTime;
result.pgf = @(z, Az) geoxgeo1_pgf(z, Az, s, boundaryEmptyProb, convention);
result.analyzer = 'dqsys_geoxgeo1';
end

function Pz = geoxgeo1_pgf(z, Az, s, p0, convention)
if ~isnumeric(z) || ~isscalar(z) || ~isreal(z) || z <= 0 || z > 1
    line_error('dqsys_geoxgeo1', 'pgf argument z must lie in (0,1]');
end
if z == 1
    Pz = 1;
    return
end
denom = z - Az * (s + (1 - s) * z);
boundary = p0 * s * (z - 1) * Az / denom;
if strcmp(convention, 'LAS_DA')
    Pz = boundary;
else
    Pz = boundary * (s / z + 1 - s) + p0 * s * (1 - 1 / z);
end
end
