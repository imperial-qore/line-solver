function result = qsys_geogeo1(a, s, convention)
% QSYS_GEOGEO1 Exact closed-form analysis of a discrete-time Geo/Geo/1 queue.
%
% RESULT = QSYS_GEOGEO1(A, S) analyzes a slotted single-server queue with:
%   A - per-slot arrival probability (0 < A < S)
%   S - per-slot service completion probability (0 < S <= 1)
%
% RESULT = QSYS_GEOGEO1(A, S, CONVENTION) selects the observation epoch,
% 'LAS_DA' (default) or 'EAS'. Time advances in slots of unit length. In each
% slot an arrival occurs with probability A and, if the server is engaged, a
% service completion occurs with probability S, independently of everything
% else. The system content at slot boundaries is a discrete birth-death chain
% whose stationary distribution is geometric.
%
% The two conventions are not two systems. They are one system, Daduna's
% LA-rule (events at the end of their slot) with the D/A-rule (departure
% resolved before arrival), observed at two instants. With
% X(t+1) = X(t) - D(t) + A(t), 'LAS_DA' is the law of X, taken after both
% events, and 'EAS' is the law of Y(t) = X(t) - D(t), taken after the departure
% and before the arrival. They are one departure apart, so the mean contents
% differ by exactly A and, by Little's law, the sojourn times by one slot. The
% queueing delay is the same under both.
%
% With rho = A/S and r = A(1-S)/(S(1-A)):
%   LAS_DA: p_0 = 1-rho,  p_n = (1-rho) rho/(1-A) r^(n-1),  n >= 1
%           E[N] = A(1-A)/(S-A),  E[T] = (1-A)/(S-A)
%   EAS:    p_n = (1-r) r^n, n >= 0
%           E[N] = A(1-S)/(S-A),  E[T] = (1-S)/(S-A)
% Both give U = rho, X = A and E[W] = A(1-S)/(S(S-A)).
%
% Returns a struct with fields:
%   convention       - Observation epoch used, 'LAS_DA' or 'EAS'
%   arrivalProb      - Per-slot arrival probability A
%   serviceProb      - Per-slot service completion probability S
%   utilization      - Fraction of slots with the server serving, A/S
%   throughput       - Departures per slot, A
%   emptyProb        - Probability the system is empty at the epoch
%   ratio            - Geometric decay ratio of the queue-length tail, r
%   meanQueueLength  - Mean number of jobs in the system
%   meanWaitingQueue - Mean number of jobs waiting, i.e. not in service
%   meanSojournTime  - Mean sojourn time in slots
%   meanWaitingTime  - Mean waiting time in slots (epoch independent)
%   meanServiceTime  - Mean service time, 1/S under LAS_DA, (1-S)/S under EAS
%   pmf              - Function handle p(n) for the stationary queue length
%   analyzer         - Identifier string
%
% Examples:
%   r = qsys_geogeo1(0.2, 0.5);
%   r.meanQueueLength    % 0.5333
%   r.meanSojournTime    % 2.6667
%   r.pmf(0)             % 0.6
%
% Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
% Springer 2001, corollary 2.7 (the LAS_DA branch reproduces it term for term
% with b = A, p = S, c = 1-A, q = 1-S).
%
% See also QSYS_GEOXGEO1, QSYS_MM1
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(convention)
    convention = 'LAS_DA';
end
if ~ischar(convention) && ~isstring(convention)
    line_error(mfilename, 'convention must be the string ''LAS_DA'' or ''EAS''');
end
convention = upper(char(convention));
if ~strcmp(convention, 'LAS_DA') && ~strcmp(convention, 'EAS')
    line_error(mfilename, 'convention must be ''LAS_DA'' or ''EAS''');
end
if ~isnumeric(a) || ~isscalar(a) || ~isreal(a) || a <= 0 || a > 1
    line_error(mfilename, 'a must be a real scalar in (0,1]');
end
if ~isnumeric(s) || ~isscalar(s) || ~isreal(s) || s <= 0 || s > 1
    line_error(mfilename, 's must be a real scalar in (0,1]');
end
if a >= s
    line_error(mfilename, 'load a/s must be strictly less than 1');
end

rho = a / s;
ratio = a * (1 - s) / (s * (1 - a));

% The queueing delay does not depend on the observation epoch; only the
% accounting of the slot in which service takes place does.
meanWaitingTime = a * (1 - s) / (s * (s - a));
meanWaitingQueue = a * meanWaitingTime;

if strcmp(convention, 'LAS_DA')
    emptyProb = 1 - rho;
    meanQueueLength = a * (1 - a) / (s - a);
    meanSojournTime = (1 - a) / (s - a);
    meanServiceTime = 1 / s;
    pmf = @(n) las_da_pmf(n, emptyProb, rho, a, ratio);
else
    emptyProb = 1 - ratio;
    meanQueueLength = a * (1 - s) / (s - a);
    meanSojournTime = (1 - s) / (s - a);
    meanServiceTime = (1 - s) / s;
    pmf = @(n) eas_pmf(n, ratio);
end

result = struct();
result.convention = convention;
result.arrivalProb = a;
result.serviceProb = s;
result.utilization = rho;
result.throughput = a;
result.emptyProb = emptyProb;
result.ratio = ratio;
result.meanQueueLength = meanQueueLength;
result.meanWaitingQueue = meanWaitingQueue;
result.meanSojournTime = meanSojournTime;
result.meanWaitingTime = meanWaitingTime;
result.meanServiceTime = meanServiceTime;
result.pmf = pmf;
result.analyzer = 'qsys_geogeo1';
end

function pr = las_da_pmf(n, emptyProb, rho, a, ratio)
if any(n < 0) || any(floor(n) ~= n)
    line_error('qsys_geogeo1', 'queue length must be a non-negative integer');
end
pr = zeros(size(n));
pr(n == 0) = emptyProb;
k = n(n >= 1);
pr(n >= 1) = emptyProb * (rho / (1 - a)) * ratio .^ (k - 1);
end

function pr = eas_pmf(n, ratio)
if any(n < 0) || any(floor(n) ~= n)
    line_error('qsys_geogeo1', 'queue length must be a non-negative integer');
end
pr = (1 - ratio) * ratio .^ n;
end
