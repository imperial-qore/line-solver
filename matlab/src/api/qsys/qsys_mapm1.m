function result = qsys_mapm1(D0, D1, mu, varargin)
% QSYS_MAPM1 Analyzes a MAP/M/1 queue using Q-MAM.
%
% RESULT = QSYS_MAPM1(D0, D1, MU) analyzes a MAP/M/1 queue with:
%   D0 - MAP hidden transition matrix (n x n)
%   D1 - MAP arrival transition matrix (n x n)
%   MU - Exponential service rate
%
% This is a convenience wrapper for QSYS_MAPMC with c=1.
%
% RESULT = QSYS_MAPM1(..., 'maxNumComp', N) sets max queue length probs (default 500)
%
% Returns a struct with fields:
%   meanQueueLength   - Mean number of customers in system
%   meanWaitingTime   - Mean waiting time in queue
%   meanSojournTime   - Mean sojourn time (waiting + service)
%   utilization       - Server utilization
%   queueLengthDist   - Queue length distribution P(Q=n)
%   waitingTimePH     - Struct with alpha and T for waiting time PH
%   analyzer          - Name of analyzer used
%
% See also QSYS_MAPMC, qsys_mapph1

result = qsys_mapmc(D0, D1, mu, 1, varargin{:});
result.analyzer = 'Q-MAM:MAP/M/1';

end
