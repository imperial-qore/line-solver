function result = dqsys_bernoulli1(b, p, L)
% DQSYS_BERNOULLI1 Exact analysis of a state dependent Bernoulli server.
%
% RESULT = DQSYS_BERNOULLI1(B, P) analyzes the discrete-time single-server
% queue of Daduna (2001), chapter 2, with an unbounded buffer:
%   B - per-slot arrival probability, scalar b, or vector with B(n+1)=b(n)
%   P - per-slot service completion probability, scalar p, or vector with
%       P(n)=p(n) for n = 1,2,...
%
% RESULT = DQSYS_BERNOULLI1(B, P, L) caps the buffer at L jobs. Arrivals in a
% slot that finds L jobs present are lost, which is the loss system of
% corollary 2.8. An unbounded buffer requires scalar B and P with B < P.
%
% Time advances in slots. In the slot starting at t with n jobs present, the
% job in service departs with probability p(n) and an arrival occurs with
% probability b(n), independently; both are recorded at the end of the slot
% with the departure resolved first (Daduna's LA rule and D/A rule). The queue
% length at slot boundaries is a discrete birth-death chain with
%
%   pi(n) = [prod_{m=0}^{n-1} b(m) / prod_{m=0}^{n} c(m)]
%         * [prod_{m=1}^{n-1} q(m) / prod_{m=1}^{n} p(m)] / H,
%
% c = 1-b and q = 1-p, which is theorem 2.3 (and corollary 2.8 once b(n)=0
% above the capacity). For constant b and p it collapses to the Geo/Geo/1 law
% of DQSYS_GEOGEO1 under the LAS_DA convention.
%
% The distribution seen by an arriving customer, with himself not counted, is
% theorem 2.11,
%
%   pi_1(n) = [prod_{m=0}^{n} b(m) / prod_{m=0}^{n+1} c(m)]
%           * [prod_{m=1}^{n} q(m) / prod_{m=1}^{n} p(m)] / H_1,
%
% and is returned in ARRIVALPMF. It is not the time-stationary law: discrete
% time has no PASTA analogue, and the two laws differ even when the arrival
% stream is a state independent Bernoulli process. For that state independent
% case pi_1 is exactly the EAS-convention queue length law of DQSYS_GEOGEO1,
% geometric with ratio r = b(1-p)/(p(1-b)).
%
% Returns a struct with fields:
%   capacity         - Buffer capacity, Inf when unbounded
%   arrivalProb      - Offered per-slot arrival probability, b(n)
%   serviceProb      - Per-slot service completion probability, p(n)
%   pmf              - Time-stationary queue length law, index n+1 holds pi(n)
%   arrivalPmf       - Arrival queue length law of theorem 2.11
%   emptyProb        - pi(0)
%   utilization      - Fraction of slots with the server busy, 1 - pi(0)
%   throughput       - Carried departures per slot
%   lossProb         - Fraction of offered arrivals lost, 0 when unbounded
%   meanQueueLength  - Mean number of jobs in the system
%   meanWaitingQueue - Mean number of jobs waiting, i.e. not in service
%   meanSojournTime  - Mean sojourn time in slots, by Little's law
%   meanWaitingTime  - Mean waiting time in slots
%   normConst        - Normalizing constant H of theorem 2.3
%   analyzer         - Identifier string
%
% Examples:
%   r = dqsys_bernoulli1(0.2, 0.5);          % Geo/Geo/1
%   r.meanQueueLength                        % 0.5333
%   r = dqsys_bernoulli1(0.2, 0.5, 4);       % Geo/Geo/1/4 loss system
%   r.lossProb
%   r = dqsys_bernoulli1(0.6, 0.3*min(1:20,3), 20);   % example 2.10, s = 3
%
% Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
% Springer 2001, theorem 2.3, corollaries 2.7 and 2.8, example 2.10 and
% theorem 2.11.
%
% See also DQSYS_GEOGEO1, DQSYS_GEOXGEO1, DPFQN_NC
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(L)
    L = Inf;
end
if ~isnumeric(L) || ~isscalar(L) || ~isreal(L) || L < 1 || (isfinite(L) && floor(L) ~= L)
    line_error(mfilename, 'L must be a positive integer or Inf');
end
if ~isnumeric(b) || ~isreal(b) || isempty(b) || any(b(:) < 0) || any(b(:) > 1)
    line_error(mfilename, 'arrival probabilities must be real and in [0,1]');
end
if ~isnumeric(p) || ~isreal(p) || isempty(p) || any(p(:) <= 0) || any(p(:) > 1)
    line_error(mfilename, 'service probabilities must be real and in (0,1]');
end

if isinf(L)
    if ~isscalar(b) || ~isscalar(p)
        line_error(mfilename, 'an unbounded buffer requires scalar arrival and service probabilities; pass a capacity L to analyze a state dependent server');
    end
    if b >= p
        line_error(mfilename, 'load b/p must be strictly less than 1 on an unbounded buffer');
    end
    result = bernoulli1_infinite(b, p);
    return
end

% Finite capacity: expand b and p over the reachable states 0..L. The offered
% stream is kept separately from the admitted one, so the loss probability
% stays meaningful for a state dependent offered stream.
boff = expand_state(b, L+1, 'arrival');       % boff(n+1) = b(n), n = 0..L
pv = expand_service(p, L);                    % pv(n) = p(n), n = 1..L
badm = boff;
badm(L+1) = 0;                                % an arrival finding L jobs is lost
if any(badm(1:L) >= 1)
    % c(n)=0 makes the weight of theorem 2.3 diverge at n; p(n)=1 is fine and
    % truncates the chain instead, which is the deterministic service of
    % example 2.9.
    line_error(mfilename, 'arrival probabilities below the capacity must be strictly less than one');
end

lw = zeros(1, L+1);                           % log unnormalized pi, theorem 2.3
acc = -log(1 - badm(1));
lw(1) = acc;
for n = 1:L
    acc = acc + log(badm(n)) - log(1 - badm(n+1)) - log(pv(n));
    if n >= 2
        acc = acc + log(1 - pv(n-1));
    end
    lw(n+1) = acc;
end
w = exp(lw - max(lw));
H = sum(w);
pmf = w / H;

la = zeros(1, L);                             % log unnormalized pi_1, theorem 2.11
if L >= 1
    acc = log(badm(1)) - log(1 - badm(1)) - log(1 - badm(2));
    la(1) = acc;
    for n = 1:L-1
        acc = acc + log(badm(n+1)) - log(1 - badm(n+2)) + log(1 - pv(n)) - log(pv(n));
        la(n+1) = acc;
    end
end
if isempty(la) || all(~isfinite(la))
    arrivalPmf = zeros(1, 0);
else
    wa = exp(la - max(la(isfinite(la))));
    wa(~isfinite(la)) = 0;
    arrivalPmf = wa / sum(wa);
end

nvec = 0:L;
meanQueueLength = sum(pmf .* nvec);
utilization = 1 - pmf(1);
throughput = sum(pmf(2:end) .* pv(1:L));
offered = sum(pmf .* boff);
if offered > 0
    lossProb = pmf(L+1) * boff(L+1) / offered;
else
    lossProb = 0;
end
meanWaitingQueue = meanQueueLength - utilization;

result = struct();
result.capacity = L;
result.arrivalProb = boff;
result.serviceProb = pv;
result.pmf = pmf;
result.arrivalPmf = arrivalPmf;
result.emptyProb = pmf(1);
result.utilization = utilization;
result.throughput = throughput;
result.lossProb = lossProb;
result.meanQueueLength = meanQueueLength;
result.meanWaitingQueue = meanWaitingQueue;
if throughput > 0
    result.meanSojournTime = meanQueueLength / throughput;
    result.meanWaitingTime = meanWaitingQueue / throughput;
else
    result.meanSojournTime = 0;
    result.meanWaitingTime = 0;
end
result.normConst = H * exp(max(lw));
result.analyzer = 'dqsys_bernoulli1';
end

% ==========================================================================
function res = bernoulli1_infinite(b, p)
% Corollary 2.7 in closed form, with the arrival law of theorem 2.11.
g = dqsys_geogeo1(b, p, 'LAS_DA');
r = g.ratio;
res = struct();
res.capacity = Inf;
res.arrivalProb = b;
res.serviceProb = p;
res.pmf = @(n) g.pmf(n);
res.arrivalPmf = @(n) arrival_geom(n, r);
res.emptyProb = g.emptyProb;
res.utilization = g.utilization;
res.throughput = g.throughput;
res.lossProb = 0;
res.meanQueueLength = g.meanQueueLength;
res.meanWaitingQueue = g.meanWaitingQueue;
res.meanSojournTime = g.meanSojournTime;
res.meanWaitingTime = g.meanWaitingTime;
res.normConst = 1 / g.emptyProb;
res.analyzer = 'dqsys_bernoulli1';
end

% ==========================================================================
function pr = arrival_geom(n, r)
if any(n < 0) || any(floor(n) ~= n)
    line_error('dqsys_bernoulli1', 'queue length must be a non-negative integer');
end
pr = (1 - r) * r .^ n;
end

% ==========================================================================
function v = expand_state(x, len, what)
% Expand a scalar or a vector indexed by n = 0..len-1 to full length.
if isscalar(x)
    v = repmat(x, 1, len);
elseif numel(x) == len
    v = x(:)';
else
    line_error('dqsys_bernoulli1', 'the %s probability vector must have %d entries, one per state 0..%d', what, len, len-1);
end
end

% ==========================================================================
function v = expand_service(x, L)
% Expand the service probabilities to p(1..L). A vector of length L+1 is
% accepted with its first entry, which would be p(0), ignored.
if isscalar(x)
    v = repmat(x, 1, L);
elseif numel(x) == L
    v = x(:)';
elseif numel(x) == L + 1
    x = x(:)';
    v = x(2:end);
else
    line_error('dqsys_bernoulli1', 'the service probability vector must have %d or %d entries', L, L+1);
end
end
