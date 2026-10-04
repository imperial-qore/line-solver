function result = qsys_mmk_qed(lambda, mu, s)
% QSYS_MMK_QED Halfin-Whitt QED approximation for the M/M/s queue.
%
% RESULT = QSYS_MMK_QED(LAMBDA, MU, S) returns the quality-and-efficiency-driven
% (QED) approximation for the M/M/S queue with arrival rate LAMBDA and service
% rate MU per server.
%
% THE REGIME. Let S grow with the offered load a = LAMBDA/MU so that the SERVER
% SLACK stays of order sqrt(S), i.e.
%
%   beta = (1 - rho) sqrt(S) = (S - a)/sqrt(S)
%
% is held fixed. Halfin and Whitt proved that the delay probability then has the
% non-degenerate limit
%
%   alpha(beta) = [ 1 + beta Phi(beta)/phi(beta) ]^(-1)
%
% with phi and Phi the standard normal density and cdf. That is the whole point
% of the regime: servers are busy a fraction 1 - beta/sqrt(S) of the time, so
% efficiency tends to 1, and yet the delay probability tends to a constant
% strictly between 0 and 1, so quality does not collapse. Neither the
% underloaded regime (alpha -> 0) nor the overloaded one (alpha -> 1) has that
% property.
%
% WHY USE IT WHEN M/M/s IS EXACT. Erlang C needs a sum of S terms of the form
% a^j/j!, which overflows in double precision well before the thousands of
% servers a large contact centre or a datacentre thread pool has; alpha(beta) is
% three transcendental calls at any S. The approximation is also the object that
% the staffing rule inverts, see QSYS_MMK_QED_STAFFING.
%
% Returns a struct with fields:
%   offeredLoad      - a = lambda/mu, in erlangs
%   trafficIntensity - rho = a/s
%   beta             - the QED server-slack parameter (s-a)/sqrt(s)
%   probDelay        - alpha(beta), the probability an arrival waits
%   meanWaitDelayed  - E[W | W > 0] = 1/(s*mu - lambda), exact for M/M/s
%   meanWait         - E[W] = alpha(beta)/(s*mu - lambda)
%   meanQueueLength  - E[Q] = lambda E[W], customers waiting
%   meanNumber       - E[N] = a + E[Q]
%   utilization      - rho
%
% ACCURACY. The error is O(1/sqrt(s)) at fixed beta: at s = 100 and beta = 0.5
% the exact Erlang C is 0.5065 against 0.5045 here, and at s = 100000 it is
% 0.50461 against 0.50454.
%
% An overloaded model, beta <= 0, has no QED limit; probDelay is then reported
% as 1 and the waiting-time fields as Inf, which is what the M/M/s queue does.
%
% Example:
%   res = qsys_mmk_qed(990, 1, 1000);   % 1000 servers, 99% loaded
%   res.probDelay
%
% Reference: S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with
% many exponential servers. Operations Research 29(3), 567-588.
%
% See also QSYS_MMK, QSYS_MMK_QED_STAFFING, QSYS_ERLANGA.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if lambda <= 0
    line_error(mfilename, 'The arrival rate lambda must be positive.');
end
if mu <= 0
    line_error(mfilename, 'The service rate mu must be positive.');
end
s = round(s);
if s < 1
    line_error(mfilename, 'The number of servers s must be at least 1.');
end

a = lambda / mu;
rho = a / s;
beta = (s - a) / sqrt(s);

result.offeredLoad = a;
result.trafficIntensity = rho;
result.beta = beta;
result.utilization = rho;

if beta <= 0
    % Not a QED model: every arrival is delayed and the queue has no
    % steady state, exactly as in the M/M/s queue at rho >= 1.
    result.probDelay = 1;
    result.meanWaitDelayed = Inf;
    result.meanWait = Inf;
    result.meanQueueLength = Inf;
    result.meanNumber = Inf;
    return
end

result.probDelay = qsys_mmk_qed_alpha(beta);
result.meanWaitDelayed = 1 / (s*mu - lambda);
result.meanWait = result.probDelay * result.meanWaitDelayed;
result.meanQueueLength = lambda * result.meanWait;
result.meanNumber = a + result.meanQueueLength;
end
