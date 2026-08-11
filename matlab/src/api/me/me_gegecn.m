function [p, L, U, PB, Lq] = me_gegecn(lambda, Ca, mu, Cs, c, K, N)
%ME_GEGECN Censored GE/GE/c/K;N queue by entropy maximisation
%
% Maximum Entropy solution of the single-class censored FCFS GE/GE/c/K;N
% queue of Kouvatsos (1994), Section 4.1. The queue holds at most N jobs
% and never fewer than K; arrivals finding N jobs present are turned away
% and departures are not allowed from state K. For a queue embedded in an
% OPEN network K is always 0; positive K arises in closed networks, where
% it records the minimum occupancy forced by the remaining stations being
% full.
%
% The ME state probabilities subject to normalisation, the marginal
% utilisations u(l), the mean queue length excluding J jobs and the
% full-buffer probability coincide with the global balance solution
%
%   p(n) = p(K) * G_n * x^h(n) * y^f(n),   n = K+1,...,N          (4.2)
%
% with G_n = prod_{l=K+1}^{m(n)} g(l), J = max(c,K+1),
% h(n) = max(0,n-J), f(n) = max(0,n-N+1), m(n) = max{K+1,min(c,n)}, and
% the Lagrangian coefficients g(l), x, y given in closed form by raw
% system data. The coefficients are invariant to N and K, so letting
% K -> 0 and N -> Inf recovers the stable GE/GE/c solution used by
% ME_OQN.
%
% INPUTS:
%   lambda - Arrival rate offered to the queue (attempts, including the
%            arrivals that are turned away)
%   Ca     - Squared coefficient of variation of the interarrival times
%            (Ca >= 1: the GE distribution is defined for scv >= 1)
%   mu     - Service rate of one server
%   Cs     - Squared coefficient of variation of the service times
%            (Cs >= 1)
%   c      - Number of servers (c >= 1, finite)
%   K      - Minimum number of jobs in the queue (K >= 0)
%   N      - Buffer capacity in jobs, in service included (N > K, finite)
%
% OUTPUTS:
%   p  - Queue length distribution, p(idx) = Pr{n = K+idx-1}, idx = 1..N-K+1
%   L  - Mean number of jobs in the queue, sum_n n*p(n)
%   U  - Utilization, E[min(n,c)]/c (mean fraction of busy servers)
%   PB - Probability that an arrival of the queue's own aggregate stream
%        finds the queue full, eq. (4.3)
%   Lq - Mean number of jobs waiting, L - E[min(n,c)]
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994, Section 4.1,
%   equations (4.1)-(4.3).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isfinite(N)
    line_error(mfilename, 'me_gegecn requires a finite buffer capacity N; use me_oqn for infinite capacity.');
end
if ~isfinite(c) || c < 1
    line_error(mfilename, 'me_gegecn requires a finite number of servers c >= 1.');
end
if N <= K
    line_error(mfilename, 'me_gegecn requires N > K.');
end
if Ca < 1 - 1e-12 || Cs < 1 - 1e-12
    line_error(mfilename, 'me_gegecn requires Ca >= 1 and Cs >= 1: the GE distribution is not defined for scv < 1.');
end
if mu <= 0
    line_error(mfilename, 'me_gegecn requires a positive service rate.');
end

c = round(c);
K = round(K);
N = round(N);

tau = 2 / (Ca + 1);
sigma = 2 / (Cs + 1);
rho = lambda / (c * mu);

J = max(c, K + 1);
den1 = sigma * (1 - tau) + tau;   % sigma(1-tau)+tau
den2 = tau * rho * (1 - sigma) + sigma; % tau*rho(1-sigma)+sigma

% Lagrangian coefficients g(l), l = K+1,...,J
g = ones(J, 1);
if K < c - 1
    g(K+1) = tau * c * rho / ((K + 1) * den1);
elseif K == c - 1
    g(K+1) = tau * sigma * rho / den2;
else % K >= c
    g(K+1) = (den1 / den2) * tau * rho;
end
for l = (K+2):J
    if l < J
        g(l) = (tau * c * rho + (l - 1) * sigma * (1 - tau)) / (l * den1);
    else
        g(l) = sigma * (tau * c * rho + (J - 1) * sigma * (1 - tau)) / (J * den2);
    end
end

x = (tau * rho + sigma * (1 - tau)) / den2;
y = 1 / (1 - (1 - sigma) * x);

% Unnormalized log-probabilities. Working in log space keeps x^(N-J) from
% overflowing on a saturated queue with a large buffer, and makes the
% rho=1 case (x=1) fall out of the same expression instead of needing the
% separate p(K) branch of (4.2).
n = K:N;
cumlogg = cumsum(log(g(K+1:J)));
logp = zeros(1, numel(n));
for idx = 1:numel(n)
    nn = n(idx);
    if nn > K
        m = max(K + 1, min(c, nn));
        logp(idx) = cumlogg(m - K);
    end
end
logp = logp + max(0, n - J) * log(x) + max(0, n - N + 1) * log(y);
logp = logp - max(logp);
p = exp(logp);
p = p / sum(p);

busy = min(n, c);
L = sum(n .* p);
Ebusy = sum(busy .* p);
U = Ebusy / c;
Lq = L - Ebusy;
PB = me_gegecn_pb(p, K, N, c, Cs, Ca);
end
